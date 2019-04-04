from __future__ import print_function
from parameters import *
import numpy as np

# Define initial values
x = np.array([1.0])
x_n = np.array([1.0])
y = np.array([1.0, 0.0])
y_n = np.array([1.0, 0.0])

# Define function solving fluid subproblem
def fluid(x_n, y_n, y, save):
    x_N = np.array([x_n[0]])
    for m in range(M):
        A = M/dt + gamma/2.0
        a = M/dt*x_n[0] - gamma/2.0*x_n[0] + \
            beta/2.0*((M - m - 1.0)/M*y_n[0] + (m + 1.0)/M*y[0])* \
                     ((M - m - 1.0)/M*y_n[0] + (m + 1.0)/M*y[0]) +\
            beta/2.0*((M - m)/M*y_n[0] + m/M*y[0])* \
                     ((M - m)/M*y_n[0] + m/M*y[0])
        x = np.array([a/A])
        if (save == True):
            u_f_txt = open('u_f.txt', 'a')
            u_f_txt.write(str((M*i + m + 1)*dt/M) + ' ' + str(x[0]) + '\r\n')
            u_f_txt.close()
            print('Fluid displacement: ', str(x[0]))
        x_n[0] = x[0]
    x_n[0] = x_N[0]
    return x

# Define function solving solid subproblem
def solid(y_n, x_n, x, save):
    y_N = np.array([y_n[0], y_n[1]])
    for k in range(K):
        B = np.array([[K/dt, mu/2.0] ,\
                      [-1.0/2.0, K/dt]])
        b = np.array([K/dt*y_n[0] - mu/2.0*y_n[1] + \
                      alpha/2.0*((K - k - 1.0)/K*x_n[0] + (k + 1.0)/K*x[0]) + \
                      alpha/2.0*((K - k)/K*x_n[0] + k/K*x[0]), 
                      K/dt*y_n[1] + 1.0/2.0*y_n[0]])
        y = np.linalg.solve(B, b)
        if (save == True):
            v_s_txt = open('v_s.txt', 'a')
            v_s_txt.write(str((K*i + k + 1)*dt/K) + ' ' + str(y[0]) + '\r\n')
            v_s_txt.close()
            print('Solid velocity: ', str(y[0]))
            u_s_txt = open('u_s.txt', 'a')
            u_s_txt.write(str((K*i + k + 1)*dt/K) + ' ' + str(y[1]) + '\r\n')
            u_s_txt.close()
            print('Solid displacement: ', str(y[1]))
        y_n[0] = y[0]
        y_n[1] = y[1]
    y_n[0] = y_N[0]
    y_n[1] = y_N[1]
    return y

# Define function for shooting method
def F(s, x_n, y_n):
    y0 = np.array([s, 0.0])
    x = fluid(x_n, y_n, y0, False)
    y = solid(y_n, x_n, x, False)
    return s - y[0]

# Save and print initial
u_f_txt = open('u_f.txt', 'a')
u_f_txt.write(str(0.0) + ' ' + str(x[0]) + '\r\n')
u_f_txt.close()
print('Fluid displacement: ', str(x[0]))
v_s_txt = open('v_s.txt', 'a')
v_s_txt.write(str(0.0) + ' ' + str(y[0]) + '\r\n')
v_s_txt.close()
print('Solid velocity: ', str(y[0]))
u_s_txt = open('u_s.txt', 'a')
u_s_txt.write(str(0.0) + ' ' + str(y[1]) + '\r\n')
u_s_txt.close()
print('Solid displacement: ', str(y[1]))

# Define time loop
for i in range(N):
    
    print('Number of time step: ', str(i + 1))
    
    # Perform decoupling iterations
    s = y[0]
    x_n[0] = x[0]
    y_n[0] = y[0]
    y_n[1] = y[1]
    stop = False
    l = 0
    while (stop == False):
        
        # Perform shooting
        f = F(s, x_n, y_n)
        f_eps = F(s + eps, x_n, y_n)
        f_der = (f_eps - f)/eps
        d = -f/f_der
        s += d
        
        # Check stop conditions
        l += 1
        if (abs(f) < tol):
            stop = True
            print('Method converged successfully after ', str(l), ' iterations.')
        if (l == maxit):
            stop = True
            print('The maximal number of iterations was reached.')
            print('Value of function F: ', str(f))
            
    # Save solution and advance in time
    y0 = np.array([s, 0.0])
    x = fluid(x_n, y_n, y0, True)
    y = solid(y_n, x_n, x, True)
