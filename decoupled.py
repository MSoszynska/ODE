from __future__ import print_function
from parameters import *
import numpy as np

# Define initial values
x = np.array([1.0])
x_n = np.array([1.0])
x_new = np.array([1.0])
y = np.array([1.0, 0.0])
y_n = np.array([1.0, 0.0])

# Define time loop
for i in range(N):
    
    print('Number of time step: ', str(i + 1))
    
    # Save and print values
    u_f_txt = open('u_f.txt', 'a')
    u_f_txt.write(str(i*dt) + ' ' + str(x[0]) + '\r\n')
    u_f_txt.close()
    print('Fluid displacement: ', str(x[0]))
    v_s_txt = open('v_s.txt', 'a')
    v_s_txt.write(str(i*dt) + ' ' + str(y[0]) + '\r\n')
    v_s_txt.close()
    print('Solid velocity: ', str(y[0]))
    u_s_txt = open('u_s.txt', 'a')
    u_s_txt.write(str(i*dt) + ' ' + str(y[1]) + '\r\n')
    u_s_txt.close()
    print('Solid displacement: ', str(y[1]))
    
    # Perform decoupling iterations
    x_n[0] = x[0]
    y_n[0] = y[0]
    y_n[1] = y[1]
    stop = False
    k = 0
    while (stop == False):
        
        # Fluid step
        A = 1.0/dt + gamma/2.0
        a = 1.0/dt*x_n[0] - gamma/2.0*x_n[0] + \
            beta/2.0*y[0]*y[0] + beta/2.0*y_n[0]*y_n[0]
        x_new[0] = a/A
        
        # Solid step
        B = np.array([[1.0/dt, mu/2.0] ,\
                      [-1.0/2.0, 1.0/dt]])
        b = np.array([1.0/dt*y_n[0] - mu/2.0*y_n[1] + alpha/2*x_new[0] + alpha/2.0*x_n[0], 
                      1.0/dt*y_n[1] + 1.0/2.0*y_n[0]])
        y_new = np.linalg.solve(B, b)
        
        # Check stop conditions
        k += 1
        if (abs(x[0] - x_new[0]) < tol) and (np.linalg.norm(y - y_new) < tol):
            stop = True
            print('Method converged successfully after ', str(k), ' iterations.')
        if (k == maxit):
            stop = True
            print('The maximal number of iterations was reached.')
            print('Norm of fluid solution on the interface: ', str(abs(x[0] - x_new[0])))
            print('Norm of solid solution on the interface', str(np.linalg.norm(y - y_new)))
            
        # Relaxation
        x[0] = tau*x_new[0] + (1.0 - tau)*x[0]
        y[0] = tau*y_new[0] + (1.0 - tau)*y[0]
        y[1] = tau*y_new[1] + (1.0 - tau)*y[1]
            
# Save and print values from the last time step
u_f_txt = open('u_f.txt', 'a')
u_f_txt.write(str(N*dt) + ' ' + str(x[0]) + '\r\n')
u_f_txt.close()
print('Fluid displacement: ', str(x[0]))
v_s_txt = open('v_s.txt', 'a')
v_s_txt.write(str(N*dt) + ' ' + str(y[0]) + '\r\n')
v_s_txt.close()
print('Solid velocity: ', str(y[0]))
u_s_txt = open('u_s.txt', 'a')
u_s_txt.write(str(N*dt) + ' ' + str(y[1]) + '\r\n')
u_s_txt.close()
print('Solid displacement: ', str(y[1]))
