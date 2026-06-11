CLC-DE: Coupled Lévy-Convolutional Differential Evolution
1. Overview
CLC-DE is an improved Differential Evolution algorithm designed for continuous numerical optimization problems. 
This implementation is written in MATLAB.
________________________________________
2. File Structure
Place the following files in the same directory:
CLC_DE.m
README.md
The main algorithm file is:
CLC_DE.m
The helper functions are included inside CLC_DE.m, including:
cauchyrnd
cauchyinv
Levy
Therefore, no additional function files are required in most cases.
________________________________________
3. Function Interface
The algorithm can be called as follows:
[gbest, gbestval, con] = CLC_DE(fhd, dim, N, Max_Gen, VRmin, VRmax, varargin)
________________________________________
4. Input Arguments
Argument	Description
fhd	Function handle of the objective function
dim	Dimension of the optimization problem
N	Initial population size
Max_Gen	Maximum number of generations
VRmin	Lower bound of decision variables; it can be a scalar or a vector of length dim
VRmax	Upper bound of decision variables; it can be a scalar or a vector of length dim
varargin	Optional additional parameters passed to the objective function
________________________________________
5. Output Arguments
Output	Description
gbest	Best solution found by the algorithm, returned as a 1 × dim row vector
gbestval	Best objective function value
con	Convergence curve recording the best fitness value during the optimization process
________________________________________
6. Objective Function Format
CLC-DE is designed for minimization problems.
The objective function should support matrix-form input:
fitness = fhd(X)
where X is a dim × N matrix. Each column of X represents one candidate solution.
The objective function should return a fitness vector:
fitness
The recommended output size is 1 × N, where each element corresponds to the objective value of one candidate solution.
For example, the Sphere function can be defined as:
function fitness = sphere_func(X)
    fitness = sum(X.^2, 1);
end
________________________________________
7. Basic Usage Example
7.1 Create a Test Function
Create a file named sphere_func.m:
function fitness = sphere_func(X)
    fitness = sum(X.^2, 1);
end
7.2 Run CLC-DE
Create a script named main.m:
clc;
clear;
close all;

% Add current directory
addpath(pwd);

% Parameter settings
dim = 30;          % Problem dimension
N = 100;           % Initial population size
Max_Gen = 1000;    % Maximum number of generations
VRmin = -100;      % Lower bound
VRmax = 100;       % Upper bound

% Run CLC-DE
[gbest, gbestval, con] = CLC_DE(@sphere_func, dim, N, Max_Gen, VRmin, VRmax);

% Display results
disp('Best solution:');
disp(gbest);

disp('Best fitness value:');
disp(gbestval);

% Plot convergence curve
figure;
semilogy(con, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Best Fitness Value');
title('Convergence Curve of CLC-DE');
grid on;
After running main.m, the best solution, the best fitness value, and the convergence curve will be obtained.
________________________________________
8. Example with Additional Objective Function Parameters
If the objective function requires additional parameters, they can be passed through varargin.
For example, create a shifted Sphere function named shifted_sphere.m:
function fitness = shifted_sphere(X, shift)
    fitness = sum((X - shift(:)).^2, 1);
end
Then call CLC-DE as follows:
clc;
clear;
close all;

dim = 30;
N = 100;
Max_Gen = 1000;
VRmin = -100;
VRmax = 100;

shift = ones(dim, 1) * 5;

[gbest, gbestval, con] = CLC_DE(@shifted_sphere, dim, N, Max_Gen, VRmin, VRmax, shift);

disp('Best fitness value:');
disp(gbestval);
________________________________________
9. Boundary Settings
If all decision variables have the same lower and upper bounds, scalar bounds can be used:
VRmin = -100;
VRmax = 100;
If different variables have different bounds, vector bounds can be used:
VRmin = [-10; -50; -100];
VRmax = [10; 50; 100];
dim = 3;
In this case, the lengths of VRmin and VRmax must be consistent with dim.
________________________________________
10. Main Algorithmic Components
CLC-DE includes the following main components:
1.	Adaptive parameter memory
 	The algorithm maintains historical memories for the scaling factor F and crossover rate CR. Successful parameters are used to update the memory and guide future parameter generation.
2.	External archive
 	An external archive is used to store previously successful solutions. This mechanism helps maintain population diversity and provides additional difference vectors for mutation.
3.	Linear population size reduction
 	The population size is gradually reduced during the evolutionary process to improve computational efficiency and strengthen exploitation in later generations.
4.	Dynamic mutation strategy allocation
 	The algorithm dynamically assigns different mutation strategies to individuals according to their improvement status and spatial distribution.
5.	Convolutional Filtering Mutation
 	A convolution-inspired smoothing operation was introduced to reduce unstable perturbations and guide selected individuals toward more promising search directions.
6.	Coupled Lévy Search
 	A coupled Lévy perturbation strategy combined with a rotation-like transformation is used to enhance the ability to escape local optima and improve non-axial exploration.
________________________________________
11. Recommended Parameter Settings
A typical parameter setting is:
dim = 30;
N = 100;
Max_Gen = 1000;
VRmin = -100;
VRmax = 100;
For low-dimensional or simple problems, N and Max_Gen can be reduced.
For high-dimensional or complex multimodal problems, a larger population size or more generations may be used.
The approximate maximum number of function evaluations is:
FEsmax = Max_Gen × N
Therefore, increasing N or Max_Gen will increase the computational cost.
________________________________________

