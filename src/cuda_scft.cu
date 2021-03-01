#include"struct.h"

#include "init_cuda.h"
#include "cuda_scft.h"
#include <errno.h>

#include <typeinfo>
#include"cuda_aid.cuh"

extern void average_value(std::vector<double*> data,GPU_INFO *gpu_info,CUFFT_INFO *cufft_info);
extern void getConc(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info);
extern double Free_um(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info);


extern void scft_main_loop(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info){

	cudaEvent_t start,stop;

	float msec;

	cudaError_t error;

	
	
	error=cudaEventCreate(&start);
	error=cudaEventCreate(&stop);
	error=cudaEventCreate(&start);
	error=cudaEventCreate(&stop);
	error=cudaEventRecord(start,0);	
		
	

	double energy=Free_um(gpu_info,cufft_info);
	
	error=cudaEventRecord(stop,0);	
	cudaEventSynchronize(stop);	
			
	error=cudaEventElapsedTime(&msec,start,stop);

	if(error!=cudaSuccess) printf("fft_test did not successfully detect run time\n");
			
	printf("time=%0.10f\n",msec);
	
	
	
	
	

}


extern double Free_um(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info){

	

	
	int iter=0;

	double *freeEnergy=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *freeOld=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *freeW=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *freeAB=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *freeS=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *freeDiff=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *freeWsurf=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *inCompMax=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *fpsum=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	double *psum=(double*)malloc(sizeof(double)*cufft_info->batch*gpu_info->GPU_N);
	
	double MaxInCompMax=10.0;
	do{

		MaxInCompMax=10.0;

		iter=iter+1;
	
		average_value(cufft_info->wa_cu,gpu_info,cufft_info); 

		average_value(cufft_info->wb_cu,gpu_info,cufft_info);

		getConc(gpu_info,cufft_info);

		dim3 grid(cufft_info->Nx,cufft_info->Ny,cufft_info->Nz);

		for(int gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){

			checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));	
			
			for(int i=0;i<cufft_info->NxNyNz;i++){
				//int batch=2;
				//printf("%g %g %g\n",cufft_info->Pha_cu[gpu_index][i+cufft_info->NxNyNz*batch],cufft_info->wa_cu[gpu_index][i+cufft_info->NxNyNz*batch],cufft_info->wb_cu[gpu_index][i+cufft_info->NxNyNz*batch]);
			}
			//printf("%g %g\n",cufft_info->hAB,cufft_info->lambda);
			phi_w_constrained<<<grid,cufft_info->batch,0,gpu_info->stream[gpu_index]>>>(cufft_info->wa_cu[gpu_index],cufft_info->wb_cu[gpu_index],cufft_info->pha_cu[gpu_index],cufft_info->phb_cu[gpu_index], cufft_info->Pha_cu[gpu_index],cufft_info->hAB,cufft_info->lambda,cufft_info->wopt);
	
			checkCudaErrors(cudaDeviceSynchronize());

		}
		
		if(iter%cufft_info->AverIt==0){

			MaxInCompMax=0.0;

			for(int i=0;i<cufft_info->batch*gpu_info->GPU_N;i++){	
				freeW[i]=0.0;
				freeAB[i]=0.0;
				freeS[i]=0.0;
				freeWsurf[i]=0.0;
				inCompMax[i]=0.0;
			}
			
			
		
			for(int gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++)
				for(int i=0;i<cufft_info->batch;i++){

					char fname[200];
					sprintf(fname,"./Density/pha_%d.dat",i+1);
					FILE *dp=fopen(fname,"w");
					fprintf(dp,"Nx=%d, Ny=%d, Nz=%d\n",cufft_info->Nx,cufft_info->Ny,cufft_info->Nz);
					fprintf(dp,"dx=%g, dy=%g, dz=%g\n",cufft_info->dx,cufft_info->dy,cufft_info->dz);
					
					for(long ijk=0;ijk<cufft_info->NxNyNz;ijk++){
						psum[i+gpu_index*cufft_info->batch]=1-cufft_info->pha_cu[gpu_index][ijk+i*cufft_info->NxNyNz]-cufft_info->phb_cu[gpu_index][ijk+i*cufft_info->NxNyNz];
						fpsum[i+gpu_index*cufft_info->batch]=fabs(psum[i+gpu_index*cufft_info->batch]);

						if(fpsum[i+gpu_index*cufft_info->batch]>inCompMax[i+gpu_index*cufft_info->batch]) 
							inCompMax[i+gpu_index*cufft_info->batch]=fpsum[i+gpu_index*cufft_info->batch];
						freeAB[i+gpu_index*cufft_info->batch]=freeAB[i+gpu_index*cufft_info->batch]+cufft_info->hAB*cufft_info->pha_cu[gpu_index][ijk+i*cufft_info->NxNyNz]*cufft_info->phb_cu[gpu_index][ijk+i*cufft_info->NxNyNz];
						freeW[i+gpu_index*cufft_info->batch]=freeW[i+gpu_index*cufft_info->batch]-(cufft_info->wa_cu[gpu_index][ijk+i*cufft_info->NxNyNz]*cufft_info->pha_cu[gpu_index][ijk+i*cufft_info->NxNyNz]+cufft_info->wb_cu[gpu_index][ijk+i*cufft_info->NxNyNz]*cufft_info->phb_cu[gpu_index][ijk+i*cufft_info->NxNyNz]);
						fprintf(dp,"%8.8e %8.8e %8.8e %8.8e\n",cufft_info->pha_cu[gpu_index][ijk+i*cufft_info->NxNyNz],cufft_info->phb_cu[gpu_index][ijk+i*cufft_info->NxNyNz],cufft_info->wa_cu[gpu_index][ijk+i*cufft_info->NxNyNz],cufft_info->wb_cu[gpu_index][ijk+i*cufft_info->NxNyNz]);
					}// end for ijk
					fclose(dp);
					freeAB[i+gpu_index*cufft_info->batch]/=cufft_info->NxNyNz;
					//printf("freeW=%0.10f\n",freeW);
					freeW[i+gpu_index*cufft_info->batch]/=cufft_info->NxNyNz;
					freeWsurf[i+gpu_index*cufft_info->batch]/=cufft_info->NxNyNz;
				
					freeS[i+gpu_index*cufft_info->batch]=-log(cufft_info->ql[gpu_index][i]);
					//printf("%d %.10f %.10f %.10f %.10f\n",i,qCab[0],qCab[1],freeS[i],-log(qCab[1]));
					freeOld[i+gpu_index*cufft_info->batch]=freeEnergy[i+gpu_index*cufft_info->batch];
					freeEnergy[i+gpu_index*cufft_info->batch]=freeAB[i+gpu_index*cufft_info->batch]+freeW[i+gpu_index*cufft_info->batch]+freeS[i+gpu_index*cufft_info->batch];
					printf("GPU %d batch %d: %5d : %.8e, %.8e, %.8e,%.8e, %.8e\n", gpu_index,i,iter, freeEnergy[i+gpu_index*cufft_info->batch],freeAB[i+gpu_index*cufft_info->batch],freeW[i+gpu_index*cufft_info->batch], freeS[i+gpu_index*cufft_info->batch],inCompMax[i+gpu_index*cufft_info->batch]);
					//printf("%g %g\n",inCompMax[0],inCompMax[i+gpu_index*cufft_info->batch]);
					
				

					if(inCompMax[i+gpu_index*cufft_info->batch]>MaxInCompMax) MaxInCompMax=inCompMax[i+gpu_index*cufft_info->batch];

				}// end for i
			}//end if cufft_info->AverIt==0
	
		

		//printf("%g %g\n",MaxInCompMax,cufft_info->ErrorinCompMax);

	}while(iter<cufft_info->MaxIT&&MaxInCompMax>cufft_info->ErrorinCompMax);//! end loop do


	free(freeEnergy);
	free(freeOld);
	free(freeW);
	free(freeAB);
	free(freeS);
	free(freeDiff);
	free(freeWsurf);
	free(inCompMax);
	free(psum);
	free(fpsum);
	
	return 0;

}





extern void getConc(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info){

	int gpu_index;

	dim3 grid(cufft_info->Nx,cufft_info->Ny,cufft_info->Nz);

	int threads=1024;

	size_t smemSize = threads * sizeof(double)*2;//(threads <= 32) ? 2 * threads * sizeof(double) : 

	average_value(cufft_info->wa_cu,gpu_info,cufft_info);	

	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){
		
		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));	
	
		qInt_init<<<grid,cufft_info->batch,0,gpu_info->stream[gpu_index]>>>(cufft_info->qInt_cu[gpu_index]);

		checkCudaErrors(cudaDeviceSynchronize());

	}

	sovDifFft(gpu_info,cufft_info,cufft_info->qa_cu,cufft_info->wa_cu,cufft_info->NsA,1);
	
	
	sovDifFft(gpu_info,cufft_info,cufft_info->qcb_cu,cufft_info->wb_cu,cufft_info->NsB,-1);


	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){
		
		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));	

		qa_to_qInt<<<grid,cufft_info->batch,0,gpu_info->stream[gpu_index]>>>(cufft_info->qInt_cu[gpu_index],cufft_info->qa_cu[gpu_index],cufft_info->NsA);

	}
	

	
	sovDifFft(gpu_info,cufft_info,cufft_info->qb_cu,cufft_info->wb_cu,cufft_info->NsB,1);

	
	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){
		
		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));	

		qa_to_qInt2<<<grid,cufft_info->batch,0,gpu_info->stream[gpu_index]>>>(cufft_info->qInt_cu[gpu_index],cufft_info->qcb_cu[gpu_index],cufft_info->NsB);
		
		checkCudaErrors(cudaDeviceSynchronize());

	}
	
	sovDifFft(gpu_info,cufft_info,cufft_info->qca_cu,cufft_info->wa_cu,cufft_info->NsA,-1);
	//for(int i=0;i<20;i++) printf("%g\n",cufft_info->qca_cu[0][i*cufft_info->NsA]);
	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){
		
		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));	

		cal_ql<<<cufft_info->batch,threads,smemSize,gpu_info->stream[gpu_index]>>>(cufft_info->ql[gpu_index],cufft_info->qb_cu[gpu_index],cufft_info->NsB,cufft_info->NxNyNz);

		checkCudaErrors(cudaDeviceSynchronize());
		checkCudaErrors(cudaGetLastError());

		for(int i=0;i<cufft_info->batch;i++){
			cufft_info->ql[gpu_index][i]/=cufft_info->NxNyNz;
			cufft_info->ffl[gpu_index][i]=cufft_info->ds0/cufft_info->ql[gpu_index][i];
		}
		
		//w_to_phi<<<grid,cufft_info->batch,0,gpu_info->stream[gpu_index]>>>(cufft_info->pha_cu[gpu_index], cufft_info->phb_cu[gpu_index],cufft_info->qa_cu[gpu_index],cufft_info->qca_cu[gpu_index],cufft_info->qb_cu[gpu_index],cufft_info->qcb_cu[gpu_index],cufft_info->NsA,cufft_info->NsB,cufft_info->ffl[gpu_index]);
		dim3 gridgo(cufft_info->NxNyNz/gpu_info->thread,cufft_info->batch);
		w_to_phi_go<<<gridgo,gpu_info->thread,0,gpu_info->stream[gpu_index]>>>(cufft_info->pha_cu[gpu_index], cufft_info->phb_cu[gpu_index],cufft_info->qa_cu[gpu_index],cufft_info->qca_cu[gpu_index],cufft_info->qb_cu[gpu_index],cufft_info->qcb_cu[gpu_index],cufft_info->NsA,cufft_info->NsB,cufft_info->ffl[gpu_index]);

		checkCudaErrors(cudaDeviceSynchronize());
		
	//printf("cal=%g\n",cufft_info->ql[gpu_index][0]);
	}

	checkCudaErrors(cudaGetLastError());
}



extern void sovDifFft(GPU_INFO *gpu_info,CUFFT_INFO *cufft_info,std::vector<double*> g,std::vector<double*> w,int ns,int sign){
	
	int ns1=ns+1;
	int Nx=cufft_info->Nx;
	int Ny=cufft_info->Ny;
	int Nz=cufft_info->Nz;
	int gpu_index;	
	int iz;
	
	dim3 grid(Nx,Ny,Nz),block(cufft_info->batch,1,1),grid1(cufft_info->Nxh1,Ny,Nz);
	
	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	
		
		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));
		
		initilize_wdz<<<grid,block,0,gpu_info->stream[gpu_index]>>>(w[gpu_index],cufft_info->wdz_cu[gpu_index],cufft_info->ds2);

		

	}
	
	

	if(sign==1){
		for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	

			checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));

			initilize_q<<<grid,block,1,gpu_info->stream[gpu_index]>>>(g[gpu_index],cufft_info->qInt_cu[gpu_index],ns1);//,gpu_info->stream[gpu_index]
			
		}

		for(iz=1;iz<=ns;iz++){
			for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	
		
				checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));
				dim3 gridgo(cufft_info->NxNyNz/gpu_info->thread,cufft_info->batch);

				initilize_in_go<<<gridgo,gpu_info->thread,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_in[gpu_index],g[gpu_index],cufft_info->wdz_cu[gpu_index],ns1,iz);
				//initilize_in<<<grid,block,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_in[gpu_index],g[gpu_index],cufft_info->wdz_cu[gpu_index],ns1,iz);
			
				checkCudaErrors(cufftExecD2Z(cufft_info->plan_forward[gpu_index],cufft_info->device_in[gpu_index],cufft_info->device_out[gpu_index]));
				
				
				//sufaceField<<<grid1,block,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_out[gpu_index],cufft_info->kxyzdz_cu[gpu_index],cufft_info->Nx);
				dim3 gridgo_sur(cufft_info->Nxh1NyNz/gpu_info->thread_sur,cufft_info->batch);

				sufaceField_go<<<gridgo_sur,gpu_info->thread_sur,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_out[gpu_index],cufft_info->kxyzdz_cu[gpu_index],cufft_info->Nxh1,cufft_info->Nx,cufft_info->Ny,cufft_info->Nz);
		
				checkCudaErrors(cufftExecZ2D(cufft_info->plan_backward[gpu_index],cufft_info->device_out[gpu_index],cufft_info->device_in[gpu_index]));
				
				in_to_g_go<<<gridgo,gpu_info->thread,0,gpu_info->stream[gpu_index]>>>(g[gpu_index],cufft_info->wdz_cu[gpu_index],cufft_info->device_in[gpu_index],ns1, iz);
				//in_to_g<<<grid,block,0,gpu_info->stream[gpu_index]>>>(g[gpu_index],cufft_info->wdz_cu[gpu_index],cufft_info->device_in[gpu_index],ns1, iz);
				checkCudaErrors(cudaDeviceSynchronize());
				checkCudaErrors(cudaGetLastError()); 	
				
			}
	

		}
		for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	
	
			checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));
			
			checkCudaErrors(cudaStreamSynchronize(gpu_info->stream[gpu_index]));

			checkCudaErrors(cudaDeviceSynchronize());
			
		}
		
	}
	else if(sign==-1){
		
	
		for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	
	
			checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));

			initilize_q_inverse<<<grid,block,0,gpu_info->stream[gpu_index]>>>(g[gpu_index],cufft_info->qInt_cu[gpu_index],ns1);//,gpu_info->stream[gpu_index]
			
			
			
		}
		
		
	
		for(iz=ns-1;iz>=0;iz--){
		
			for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	
	
				checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));
				dim3 gridgo(cufft_info->NxNyNz/gpu_info->thread,cufft_info->batch);
				//initilize_in<<<grid,block,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_in[gpu_index],g[gpu_index],cufft_info->wdz_cu[gpu_index],ns1,iz+2);
				initilize_in_go<<<gridgo,gpu_info->thread,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_in[gpu_index],g[gpu_index],cufft_info->wdz_cu[gpu_index],ns1,iz+2);
				
				checkCudaErrors(cufftExecD2Z(cufft_info->plan_forward[gpu_index],cufft_info->device_in[gpu_index],cufft_info->device_out[gpu_index]));

				//sufaceField<<<grid1,block,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_out[gpu_index],cufft_info->kxyzdz_cu[gpu_index],cufft_info->Nx);
				dim3 gridgo_sur(cufft_info->Nxh1NyNz/gpu_info->thread_sur,cufft_info->batch);

				sufaceField_go<<<gridgo_sur,gpu_info->thread_sur,0,gpu_info->stream[gpu_index]>>>(cufft_info->device_out[gpu_index],cufft_info->kxyzdz_cu[gpu_index],cufft_info->Nxh1,cufft_info->Nx,cufft_info->Ny,cufft_info->Nz);
		
				checkCudaErrors(cufftExecZ2D(cufft_info->plan_backward[gpu_index],cufft_info->device_out[gpu_index],cufft_info->device_in[gpu_index]));

				//in_to_g<<<grid,block,0,gpu_info->stream[gpu_index]>>>(g[gpu_index],cufft_info->wdz_cu[gpu_index],cufft_info->device_in[gpu_index],ns1, iz);
				in_to_g_go<<<gridgo,gpu_info->thread,0,gpu_info->stream[gpu_index]>>>(g[gpu_index],cufft_info->wdz_cu[gpu_index],cufft_info->device_in[gpu_index],ns1, iz);
				checkCudaErrors(cudaDeviceSynchronize());
				checkCudaErrors(cudaGetLastError()); 	
			}

		}
		for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	
	
			checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));
			
			checkCudaErrors(cudaStreamSynchronize(gpu_info->stream[gpu_index]));
			
			checkCudaErrors(cudaDeviceSynchronize());
		}
		
		
	}
	
	
	
	
	

}





extern void average_value(std::vector<double*> data,GPU_INFO *gpu_info,CUFFT_INFO *cufft_info){

	int gpu_index;
	
	int threads=gpu_info->thread;
	
	size_t smemSize = (threads <= 32) ? 2 * threads * sizeof(double) : threads * sizeof(double);

	std::vector<double*> sum;

	sum.resize(gpu_info->GPU_N);

		
	
	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	

		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));

		checkCudaErrors(cudaMallocManaged((void**)&(sum[gpu_index]), sizeof(double)* cufft_info->batch));
		
		reduce3<double><<< cufft_info->batch, threads, smemSize,gpu_info->stream[gpu_index] >>>(data[gpu_index], sum[gpu_index], cufft_info->NxNyNz);
	
		
	}
	
	
	dim3 block(cufft_info->Nx,cufft_info->Ny,cufft_info->Nz);
	
	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){
		
		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));	
		
		minus_average<<<block,cufft_info->batch,0,gpu_info->stream[gpu_index]>>>(data[gpu_index],sum[gpu_index]);
		
		
	}
	for(gpu_index=0;gpu_index<gpu_info->GPU_N;gpu_index++){	

		checkCudaErrors(cudaSetDevice(gpu_info->whichGPUs[gpu_index]));	

		checkCudaErrors(cudaFree(sum[gpu_index]));
		
		checkCudaErrors(cudaDeviceSynchronize());
	}	
	
}











