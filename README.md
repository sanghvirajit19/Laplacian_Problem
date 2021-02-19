# Laplacian problem 

![laplacian eq](https://user-images.githubusercontent.com/42026685/106930678-1c7fa180-6716-11eb-850b-84f8b7a2a6ed.png)

Above Laplacian equation is solved on the below 2-Dimentional case using Jacobi method.

![2D problem](https://user-images.githubusercontent.com/42026685/106932323-e3e0c780-6717-11eb-9a90-f60bed49b47e.png)

Jacobi step will be as followed,

![jacobi step](https://user-images.githubusercontent.com/42026685/106935413-99f9e080-671b-11eb-902c-6afe3f8fd719.png)

# Boundary condition

Boundary condition on the coarser grid of 32 x 32 is shown below.

![Problem with BC on coarser grid(32)](https://user-images.githubusercontent.com/42026685/106935942-420fa980-671c-11eb-9f8f-8878ccb3a2f2.png)

# Results

Same problem is solved on 3 different grids.
1. 32 x 32 grid
2. 64 x 64 grid
3. 128 x 128 grid

Approximate Solutions are shown below.

1. Solution on coarser grid (32x32)

    Simulation runtime: 15.942 s
    
![Approximate solutionon coarser grid (32) in 10000 Iterations](https://user-images.githubusercontent.com/42026685/106936231-9dda3280-671c-11eb-837b-22735de8cab3.png)

2. Solution on coarser grid (64x64)
  
    Simulation runtime: 68.652 s
    
![Approximate solutionon coarser grid (64) in 10000 Iterations](https://user-images.githubusercontent.com/42026685/106936690-2eb10e00-671d-11eb-8f92-4f9d6b09aba7.png)

3. Solution on finer grid (128x128)
    
    Simulation runtime: 297.846 s
    
![Approximate solutionon coarser grid (128) in 10000 Iterations](https://user-images.githubusercontent.com/42026685/106936768-425c7480-671d-11eb-8ce5-064947efbfac.png)

# Residual plot

![image](https://user-images.githubusercontent.com/42026685/108340240-a6307400-71d8-11eb-9332-f599d52c762e.png)

# Multigrid Method

![multigrid](https://user-images.githubusercontent.com/42026685/107789334-bddba880-6d51-11eb-9275-19ebb3d706e3.png)

