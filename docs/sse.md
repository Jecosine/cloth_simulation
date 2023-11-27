The code is vectorized by SSE in `cloth_code_sse.cpp`, and the auto-vectorized version is `cloth_code_vect_omp.cpp`. The manually vectorized sse version clearly outperforms all others, showcasing the importance of effective vectorization.

However, we see SSE version has the lowest IPC comparing to the other implementations. This might because:

- Vectorized operations, especially using wide AVX256 instructions, can saturate the memory bandwidth, especially if the data isn't already present in the cache. If the CPU is waiting for data to be fetched from memory, it'll stall, reducing the IPC.
- Vectorization might introduce longer dependency chains between instructions. If one instruction depends on the result of another, it can't be executed until the previous one completes.

So we find the total instructions of the sse version is much lower than the non-vectorized version. Vectorization processes multiple data points with a single instruction, reducing the overall instruction count. Its curve remains well below the other versions, indicating it's the most efficient in terms of instruction count. 