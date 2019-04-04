from __future__ import print_function
from parameters import *
import numpy as np

# Define initial values
x = np.array([1.0, 1.0, 0.0])

# Define time loop
for i in range(N):
    
    print('Number of time step: ', str(i + 1))
    
    # Save and print values
    u_f_txt = open('u_f.txt', 'a')
    u_f_txt.write(str(i*dt) + ' ' + str(x[0]) + '\r\n')
    u_f_txt.close()
    print('Fluid displacement: ', str(x[0]))
    v_s_txt = open('v_s.txt', 'a')
    v_s_txt.write(str(i*dt) + ' ' + str(x[1]) + '\r\n')
    v_s_txt.close()
    print('Solid velocity: ', str(x[1]))
    u_s_txt = open('u_s.txt', 'a')
    u_s_txt.write(str(i*dt) + ' ' + str(x[2]) + '\r\n')
    u_s_txt.close()
    print('Solid displacement: ', str(x[2]))
    
    # Perform Newton method
    x_n = np.array([x[0], x[1], x[2]])
    stop = False
    k = 0
    while (stop == False):
        F = np.array([x[0] - x_n[0] + 0.5*dt*gamma*x_n[0] - 0.5*dt*beta*x_n[1]*x_n[1] + \
                                      0.5*dt*gamma*x[0] - 0.5*dt*beta*x[1]*x[1], \
                      x[1] - x_n[1] + 0.5*dt*mu*x[2] - 0.5*dt*alpha*x[0] + \
                                      0.5*dt*mu*x[2] - 0.5*dt*alpha*x[0], \
                      x[2] - x_n[2] - 0.5*dt*x[1] - 0.5*dt*x[1]])
        J = np.array([[1.0 + dt*gamma, -2.0*dt*beta*x[1], 0.0], \
                      [-dt*alpha, 1.0, dt*mu], \
                      [0.0, -dt, 1.0]])
        d = np.linalg.solve(J, -1.0*F)
        x += d
        k += 1
        if (np.linalg.norm(F) < tol):
            print('Newton method converged successfully after ', str(k), ' iterations.')
            stop = True
        if (k == maxit):
            print('The maximal number of Newton iterations was reached.')
            print('Norm of F: ', str(np.linalg.norm(F)))
            
    
# Save and print values from the last time step
u_f_txt = open('u_f.txt', 'a')
u_f_txt.write(str(N*dt) + ' ' + str(x[0]) + '\r\n')
u_f_txt.close()
print('Fluid displacement: ', str(x[0]))
v_s_txt = open('v_s.txt', 'a')
v_s_txt.write(str(N*dt) + ' ' + str(x[1]) + '\r\n')
v_s_txt.close()
print('Solid velocity: ', str(x[1]))
u_s_txt = open('u_s.txt', 'a')
u_s_txt.write(str(N*dt) + ' ' + str(x[2]) + '\r\n')
u_s_txt.close()
print('Solid displacement: ', str(x[2]))
