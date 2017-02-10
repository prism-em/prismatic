#define BLOCK_SIZE 64
__global__ void timesTwo_d(double *a, int N){
int idx = blockDim.x * blockIdx.x + threadIdx.x;
if (idx < N)a[idx]*=2;
}



__host__ void timesTwo(double* a, const int& N){
double *a_d;
cudaMalloc((void**)&a_d, sizeof(double)*N);
cudaMemcpy(a_d,a,N*sizeof(double),cudaMemcpyHostToDevice);
timesTwo_d<<< (N-1)/BLOCK_SIZE+1, BLOCK_SIZE>>> (a_d, N);
cudaMemcpy(a,a_d,N*sizeof(double),cudaMemcpyDeviceToHost);
cudaFree(a_d);
cudaDeviceReset();
};

/*
#include <iostream>
#define N 2048
using namespace std;
int main(){
double *a = new double[N];
for (auto i =0; i<N; ++i)a[i]=i;
for (auto i =0; i<10; ++i)cout << "a[i] = " << a[i] <<", ";
cout << endl << "After doubling" << endl;;
timesTwo(a,N);
for (auto i =0; i<10; ++i)cout << "a[i] = " << a[i] <<", ";
cout << endl;
delete[] a;
return 0;
}
*/
