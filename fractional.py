from __future__ import print_function
from parameters import *
import numpy as np

# Define initial values
x = np.array([1.0])
x_n = np.array([1.0])
x_N = np.array([1.0])
x_new = np.array([1.0])
y = np.array([1.0, 0.0])
y_n = np.array([1.0, 0.0])
y_N = np.array([1.0, 0.0])

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
    x_n[0] = x[0]
    y_n[0] = y[0]
    y_n[1] = y[1]
    stop = False
    save = False
    l = 0
    while (stop == False):
        
        # Fluid step
        x_N[0] = x_n[0]
        for m in range(M):
            A = M/dt + gamma/2.0
            a = M/dt*x_n[0] - gamma/2.0*x_n[0] + \
                beta/2.0*((M - m - 1.0)/M*y_n[0] + (m + 1.0)/M*y[0])* \
                         ((M - m - 1.0)/M*y_n[0] + (m + 1.0)/M*y[0]) +\
                beta/2.0*((M - m)/M*y_n[0] + m/M*y[0])* \
                         ((M - m)/M*y_n[0] + m/M*y[0])
            x_new[0] = a/A
            if (save == True):
                u_f_txt = open('u_f.txt', 'a')
                u_f_txt.write(str((M*i + m + 1)*dt/M) + ' ' + str(x_new[0]) + '\r\n')
                u_f_txt.close()
                print('Fluid displacement: ', str(x_new[0]))
                stop = True
            x_n[0] = x_new[0]
        x_n[0] = x_N[0]
        
        # Solid step
        y_N[0] = y_n[0]
        y_N[1] = y_n[1]
        for k in range(K):
            B = np.array([[K/dt, mu/2.0] ,\
                          [-1.0/2.0, K/dt]])
            b = np.array([K/dt*y_n[0] - mu/2.0*y_n[1] + \
                          alpha/2.0*((K - k - 1.0)/K*x_n[0] + (k + 1.0)/K*x_new[0]) + \
                          alpha/2.0*((K - k)/K*x_n[0] + k/K*x_new[0]), 
                          K/dt*y_n[1] + 1.0/2.0*y_n[0]])
            y_new = np.linalg.solve(B, b)
            if (save == True):
                v_s_txt = open('v_s.txt', 'a')
                v_s_txt.write(str((K*i + k + 1)*dt/K) + ' ' + str(y_new[0]) + '\r\n')
                v_s_txt.close()
                print('Solid velocity: ', str(y_new[0]))
                u_s_txt = open('u_s.txt', 'a')
                u_s_txt.write(str((K*i + k + 1)*dt/K) + ' ' + str(y_new[1]) + '\r\n')
                u_s_txt.close()
                print('Solid displacement: ', str(y_new[1]))
            y_n[0] = y_new[0]
            y_n[1] = y_new[1]
        y_n[0] = y_N[0]
        y_n[1] = y_N[1]
        
        # Check stop conditions
        l += 1
        if (abs(x[0] - x_new[0]) < tol) and (np.linalg.norm(y - y_new) < tol) and (stop == False):
            save = True
            print('Method converged successfully after ', str(l), ' iterations.')
            print('In the next iteration solution will be saved.')
        if (l == maxit) and (stop == False):
            save = True
            print('The maximal number of iterations was reached.')
            print('In the next iteration solution will be saved.')
            print('Norm of fluid solution on the interface: ', str(abs(x[0] - x_new[0])))
            print('Norm of solid solution on the interface', str(np.linalg.norm(y - y_new)))
            
        # Relaxation
        x[0] = tau*x_new[0] + (1.0 - tau)*x[0]
        y[0] = tau*y_new[0] + (1.0 - tau)*y[0]
        y[1] = tau*y_new[1] + (1.0 - tau)*y[1]
