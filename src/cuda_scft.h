extern void sovDifFft(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info,std::vector<double*> g,std::vector<double*> w,int ns,int sign);
extern void scft_main_loop(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info);
