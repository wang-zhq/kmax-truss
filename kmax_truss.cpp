#include "kmax_truss.hpp"

using namespace std;

int main(int argc, char const **argv)
{
    time_t t0 = time(NULL);
    char const *dt_name;
    // 接收命令并读取文件
    if (argc == 1)
    {
        cout << "用法: kmtruss -f 图数据文件" << endl;
        exit(1);
    }
    else if (argc == 2) 
        dt_name = argv[1];
    else if (argc == 3)
        dt_name = argv[2];
    else
    {
        cout << "Too many input arguments." << endl;
        exit(1);
    }
    
    int64_t i, p;
    // 读取转化数据并统计无向边的总数
    size_t n_edge = 0;
    uint32_t *endpoints = import_and_trans(dt_name, &n_edge);

    cout << "数据转化完成，总边数为 " << n_edge << "。（程序已进行 " << time(NULL)-t0 << " 秒）" << endl << endl; 

    // 找到所有数里最大的那个，即为节点的总个数
    uint32_t n_node = endpoints[2*n_edge-1];

    cout << "顶点数统计完毕，源数据共有顶点 " << n_node << " 个。（程序已进行 " << time(NULL)-t0 << " 秒）" << endl << endl; 

    // cout << "----------------------------" << endl;
    // cout << "用顶点度对图数据进行精减" << endl;
    // 生成一个向量，标记每条边处于过滤端点出度的阶段
    int64_t *row_start = mark_rowhead(endpoints, n_node, n_edge);
    
    uint32_t *freq_count = new uint32_t [n_node];
    #pragma omp parallel for
    for (i = 0; i < n_node; i++)
    {
        freq_count[i] = row_start[i+1] - row_start[i];
    }
    
    uint32_t *kcore_mark = new uint32_t [n_edge];
    #pragma omp parallel for
    for (i = 0; i < n_edge; i++) {
        uint32_t a = freq_count[endpoints[2*i]];
        uint32_t b = freq_count[endpoints[2*i+1]];
        kcore_mark[i] = (a < b) ? a : b; 
    }

    // n_prem 为修剪后还剩下的边数
    int64_t n_prem;
    uint32_t try_floor = 8;
    uint32_t try_ceil = n_node;

    cout << "----------------------------" << endl;
    cout << "用线搜索方法寻找k-truss搜索下界" << endl;
    // 从点数入手就是为了减少数据模型，点数有理论缺陷，下面从三角形数量去进一步计算
    int64_t n_edge_sub = n_edge;
    uint32_t *edges_sub = new uint32_t [n_edge_sub*2];
    uint32_t n_node_sub;
    uint32_t n_try = pow(n_edge, 1./3);

    // 以1/3步长由高到低探底，得到 kmax的下界
    do
    {
        cal_ptfreq(endpoints, kcore_mark, freq_count, n_node, n_edge, n_try);
        n_prem = edge_prune(endpoints, kcore_mark, freq_count, n_node, n_edge, n_try);

        if (n_prem < n_try)
        {
            try_ceil = n_try;
            n_try *= 0.666;
            continue;
        }

        p = 0;
        for (i = 0; i < n_edge; i++)
        {
            if (kcore_mark[i] > n_try)
            {
                edges_sub[2*p] = endpoints[2*i];
                edges_sub[2*p+1] = endpoints[2*i+1];
                p++;
            }
        }
        n_edge_sub = p;
        n_node_sub = data_dim_zip(edges_sub, freq_count, n_node, n_edge_sub);

        n_prem = ktruss_chk(edges_sub, n_node_sub, n_edge_sub, n_try);
        cout << "检查数 " << n_try << " 完成，剩余边数为 " << n_prem << " (程序已进行 " <<  time(NULL)-t0 << " 秒)" << endl;
        
        if (n_prem > n_try)
            try_floor = n_try;
        else 
            try_ceil = n_try;

        
        n_try *= 0.666;
    }
    while (n_prem < n_try);
    cout << "kmax 搜索范围确定：[" << try_floor << " - " << try_ceil << "]。（程序已进行 " << time(NULL)-t0 << " 秒）" << endl << endl; 

    // 将得到了下界数据对原数据进行替换
    n_edge = n_prem;
    n_node = n_node_sub;
    memcpy(endpoints, edges_sub, sizeof(uint32_t)*n_edge*2);
    delete [] edges_sub;
    delete [] freq_count;
        
    cout << "----------------------------" << endl;
    cout << "用线搜索方法寻找 kmax" << endl;
    
    do
    {
        n_try = (try_floor+try_ceil)/2;

        n_prem = ktruss_chk(endpoints, n_node, n_edge, n_try);
        cout << "检查数 " << n_try << " 完成，剩余边数为 " << n_prem << " (程序已进行 " <<  time(NULL)-t0 << " 秒)" << endl;
        
        if (n_prem > n_try)
        {
            n_edge = n_prem;
            try_floor = n_try;
        }
        else
            try_ceil = n_try;
        
    }
    while ((try_ceil-try_floor) > 1);
    
    delete [] endpoints;
    delete [] kcore_mark;

    // 最后的 try_floor 即 kmax, 最后的 n_edge 即最终子图的边数。
    cout << endl;
    cout << "=======================================" << endl;
    cout << "kmax = " << try_floor+2 << ", Edges in kmax-truss = " << n_edge/2 << '.' << endl;
    cout << "=======================================" << endl << endl;
    cout << "ALL DONE. Time elapses " << time(NULL)-t0 << " sec." << endl;

    return 0;
}
