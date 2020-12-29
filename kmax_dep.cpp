#include "kmax_truss.hpp"

using namespace std;
// 获取文件大小
int64_t get_filesize(char const *filename)
{
    struct stat statbuf;
    int r = stat(filename, &statbuf);
    if (r == -1)
    {
        cout << "Get size of file is failed!" << endl;
        exit(-1);
    }
    return (statbuf.st_size);
}

// 读入数据并转化为二进制
uint32_t* import_and_trans(char const *dt_name, size_t *ntmp)
{
    size_t dt_size = get_filesize(dt_name);
    cout << "The size of data file " << dt_name << " is " << dt_size << " Bytes." << endl;

    // char *rawdata = (char *)malloc(sizeof(char)*dt_size);
    
    // 由于 mmap 不稳定，用 fread 读取数据
    // ifstream fin(dt_name, ios::binary);
    // fin.seekg(0L, ios::beg).read(rawdata, static_cast<streamsize>(dt_size));
    // fin.close();

    // int fd = open(dt_name, O_RDONLY);
    // int fin = read(fd, rawdata, dt_size);
    // close(fd);
    
    // FILE *fp = fopen(dt_name, "rb");
    // size_t fin = fread(rawdata, sizeof(char), dt_size, fp);
    // fclose(fp);
    
    // 通过内存映射实现快速读取
    int fp = open(dt_name, O_RDONLY);
    if (fp == -1)
    {
        printf("Can't open %s.\n", dt_name);
        exit(1);
    } 

    char *rawdata;
    rawdata = (char*)mmap(NULL, dt_size, PROT_READ, MAP_SHARED, fp, 0);
    if (rawdata == NULL || rawdata == (void*)-1)
    {
        cout << "Mapping Failed!" << endl;
        close(fp);
        exit(-2);
    }
    
    // 统计无向边的总数
    size_t n_edge = 0;
    #pragma omp parallel for reduction(+:n_edge)
    for (int64_t i = 0; i < dt_size; i++)
    {
        if (rawdata[i] == '\n') 
            n_edge++;
    }

    // 将文本数据转成二进制
    uint32_t *endpoints = new uint32_t [n_edge*2]();
    trans_txt(rawdata, endpoints, dt_size, n_edge);
    // free(rawdata);

    // 关闭镜像通道
    munmap(rawdata, dt_size);
    close(fp);

    // temptr = endpoints;
    *ntmp = n_edge;

    return endpoints;
}

// 将txt转成bin文件
// 程序已经改进，带权重和不带权重的都可以转化，用一个空格或制表符或逗号隔开的也可以识别
void trans_txt(char *rawdata, uint32_t *endpoints, size_t dt_size, size_t n_edge)
{
    uint64_t p;
    // 按线程分块，并让一块开始处于换行符后
    int n_proc = omp_get_max_threads();
    int64_t *proc_st = new int64_t [n_proc+1]();
    #pragma omp parallel for
    for (p = 1; p < n_proc; p++)
    {
        int64_t pr_end = (dt_size/n_proc)*p;
        while (rawdata[pr_end] != '\n')
            pr_end++;

        pr_end++;
        proc_st[p] = pr_end;
    }
    proc_st[n_proc] = dt_size+1;

    // 计算每一段的个数，每一段最后一位也有换行符，也得算上
    int64_t *n_dtpr = new int64_t [n_proc+1]();
    #pragma omp parallel for
    for (p = 1; p < n_proc; p++)
    {
        int64_t npz = 0;
        for (int64_t z = proc_st[p-1]; z < proc_st[p]; z++)
        {
            if (rawdata[z] == '\n')
                npz++;
        }
        n_dtpr[p] = npz*2;
    }
    
    // 算成累加序列
    for (p = 1; p <= n_proc; p++)
        n_dtpr[p] += n_dtpr[p-1];

    // 多线程并行将文本数据转换成数组
    #pragma omp parallel for
    for (p = 0; p < n_proc; p++)
    {
        int64_t pvec_z = n_dtpr[p];
        char loc_flag = 0;
        uint32_t trans_int = 0;
        for (int64_t z = proc_st[p]; z < proc_st[p+1]; z++)
        {
            if (rawdata[z] < '0')
            {
                loc_flag++;
                if ((loc_flag == 1) || (loc_flag == 2))
                    endpoints[pvec_z++] = trans_int;

                trans_int = 0;
                if (rawdata[z] == '\n')
                    loc_flag = 0;
            }
            else
            {
                trans_int *= 10;
                trans_int += (rawdata[z] - '0');
            }
        }
    }
    delete [] proc_st;
    delete [] n_dtpr; 
}

// 计算顶点出现频次
void cal_ptfreq(uint32_t *endpoints, uint32_t *kcore_mark, uint32_t *freq_count, uint32_t n_node, size_t n_edge, uint32_t n_filt)
{
    int64_t i, p;
    memset(freq_count, 0, sizeof(uint32_t)*n_node);
    #pragma omp parallel for
    for (i = 0; i < n_edge; i++)
    {
        uint32_t a = endpoints[2*i];
        uint32_t b = endpoints[2*i+1];
        if ((kcore_mark[i] > n_filt) && (a > b))
        {
            #pragma omp atomic
            freq_count[a]++; 
            #pragma omp atomic
            freq_count[b]++; 
        }
    }
}

// 去除点数出现次数不足的边
int64_t edge_prune(uint32_t *endpoints, uint32_t *kcore_mark, uint32_t *freq_count, uint32_t n_node, size_t n_edge, uint32_t n_filt)
{
    int64_t i;
    int loop_check = true;
    while (loop_check)
    {
        loop_check = false;
        
        #pragma omp parallel for
        for (i = 0; i < n_edge; i++)
        {
            if (kcore_mark[i] > n_filt)
            {
                uint32_t r_i = endpoints[i*2+1];
                uint32_t c_i = endpoints[i*2];
                if ((freq_count[c_i] <= n_filt) || (freq_count[r_i] <= n_filt))
                {
                    kcore_mark[i] = n_filt;
                    
                    #pragma omp atomic
                    freq_count[c_i]--;

                    loop_check = true;
                }
            }
        }
        cout << "\u25A1";
    }
    cout << endl;

    // 端点出现频率的总数是边的两倍
    int64_t tsm = 0;
    #pragma omp parallel for reduction(+:tsm)
    for (i = 0; i < n_node; i++)
        tsm += freq_count[i];
    
    return tsm/2;
}

// 寻找行首位置
int64_t* mark_rowhead(uint32_t *endpoints, uint32_t n_node, int64_t n_edge)
{
    int64_t *row_start_pt = new int64_t[n_node+1]();
    int64_t i;

    #pragma omp parallel for
    for (i = 1; i < n_edge; i++)
    {
        if (endpoints[i*2+1] > endpoints[(i-1)*2+1])
        {
            for (int64_t k = (endpoints[(i-1)*2+1]+1); k <= endpoints[i*2+1]; k++) 
                row_start_pt[k] = i;
        }
    }

    for (i = (endpoints[n_edge*2-1]+1); i <= n_node; i++)
        row_start_pt[i] = n_edge;

    return row_start_pt;
}

// 计算每条边所属三角形的数量
uint32_t count_tris(uint32_t *val_mx, uint32_t *endpoints, int64_t n_edge, uint32_t n_node, uint32_t n_try)
{
    int64_t i, p;

    // 定义标识点向量
    int64_t *row_start = mark_rowhead(endpoints, n_node, n_edge);

    int n_proc = omp_get_max_threads();
    char *vld = new char [n_proc*n_node]();
    
    bool drop_flag = false;
    // 有所有点对应的集合求交, 对填充率低的利用比较策略，对填充高的利用原子加
    #pragma omp parallel for
    for (i = 0; i < n_node; i++)
    {
        int64_t r_st = row_start[i];
        int64_t rl_n = row_start[i+1];

        uint32_t cid;

        int proc_id = omp_get_thread_num();
        char *vl_ini = vld + proc_id*n_node;

        int64_t k, z;
        for (k = r_st; k < rl_n; k++)
            *(vl_ini + endpoints[2*k]) = 1;

        uint32_t nti;
        for (k = r_st; k < rl_n; k++)
        {
            cid = endpoints[2*k];

            nti = 0;
            for (z = row_start[cid]; z < row_start[cid+1]; z++)
                nti += *(vl_ini + endpoints[2*z]);

            val_mx[k] = nti;
            if (nti < n_try) {
                drop_flag = true;
            }
        }

        for (k = r_st; k < rl_n; k++)
            *(vl_ini + endpoints[2*k]) = 0;

    }

    delete [] vld;
    delete [] row_start;

    if (drop_flag)
    {
        p = 0;
        for (i = 0; i < n_edge; i++)
        {
            if (val_mx[i] >= n_try)
            {
                endpoints[2*p] = endpoints[2*i];
                endpoints[2*p+1] = endpoints[2*i+1];
                p++;
            }
        }
        return p;
    }
    else
        return n_edge;

}

// 检查压缩全过程
int64_t ktruss_chk(uint32_t *endpoints, uint32_t n_node, size_t n_edge, uint32_t n_try)
{
    int64_t i, p;

    uint32_t *edges_cp = new uint32_t [n_edge*2];
    memcpy(edges_cp, endpoints, sizeof(uint32_t)*n_edge*2);
    uint32_t *val_mx = new uint32_t [n_edge];

    int64_t n_edge_cp;
    int64_t n_edge_rem = n_edge;
    do
    {
        n_edge_cp = n_edge_rem;
        n_edge_rem = count_tris(val_mx, edges_cp, n_edge_cp, n_node, n_try);

        cout << "\u25A0";
    } while ((n_edge_cp > n_edge_rem) && (n_edge_rem > n_try));
    cout << endl;

    if (n_edge_rem > n_try)
        memcpy(endpoints, edges_cp, sizeof(uint32_t)*n_edge_rem*2);
    
    delete [] edges_cp;
    delete [] val_mx;

    return n_edge_rem;
}

// 对数据进行维度压缩
uint32_t data_dim_zip(uint32_t *endpoints, uint32_t* freq_count, uint32_t n_node, size_t n_edge)
{
    int64_t i;

    uint32_t *rename_list = new uint32_t [n_node]();
    size_t n_node_cp = 0;
    for (i = 0; i < n_node; i++)
    {
        if (freq_count[i])
        {
            rename_list[i] = n_node_cp;
            n_node_cp++;
        }
    }

    // cout << "开始数据替换" << endl;
    #pragma omp parallel for
    for (i = 0; i < n_edge*2; i++)
    {
        endpoints[i] = rename_list[endpoints[i]];
    }

    delete [] rename_list;

    return n_node_cp;
}
