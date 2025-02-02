#include <iostream>
#include <cmath>

// If you have different constant double 
// You should change the number of "const double"

const double PI = 3.14; //here
const double Gravity = 9.81; //here

//Thıs ıs all of subjects and cases of Serway Mechanic

void motionandDynamics(); 
void newtonsLawofMotion(); 
void workEnergyandPower();  
void momentumandCollision();  
void rotationalMotion(); 
void angularMomentum();  
void oscillationsandWaveMotion(); 
void fluidMechanics(); 


void motionandDynamics(){
    int choice;
    double m, v, a, t, F, d, h, result;

    std::cout << "\n1 Definition of motion:(d = v0 * t + (1/2) * a * t^2)\n";
    std::cout << "2. Newton (F = m * a)\n";
    std::cout << "3. Kinetic Energy (KE = 1/2 * m * v^2)\n";
    std::cout << "4. Potential Energy (PE = m * g * h)\n";
    std::cout << "Please Enter your Selection: ";
    std::cin >> choice;

    switch (choice) {
     
        case 1:
        std::cout <<"Constant acceleration motion formula:"<< std :: endl;
        
            std::cout << "First Speed (m/s): "; std::cin >> v;
            std::cout << "Time (second): "; std::cin >> t;
            std::cout << "Angular Acceleration (m/s²): "; std::cin >> a;
            result = v * t + 0.5 * a * t * t;
            std::cout << "Final Motion: " << result << " m\n";
            break;
            // F = m.a second newton law !!!
        case 2:
        std::cout <<"Second Newtons Law Formula:"<< std :: endl;
            std::cout << "Mass (kg): "; std::cin >> m;
            std::cout << "Angular Acceleration (m/s²): "; std::cin >> a;
            result = m * a;
            std::cout << "Newton : " << result << " N\n";
            break;
            
        case 3:
        std::cout <<"Kinetic Energy:" << std :: endl;
            std::cout << "Speed (m/s): "; std::cin >> v;
            std::cout << "Mass (kg): "; std::cin >> m;
            result = 0.5 * m * v * v;
            std::cout << "Kinetic Energy: " << result << " J\n";
            break;
         
        case 4:
        std::cout <<"Gravitational Potential Energy.:"<< std :: endl;
            std::cout << "Mass(kg): "; std::cin >> m;
            std::cout << "Height (m): "; std::cin >> h;
            result = m * Gravity * h;
            std::cout << "Potential Energy is : " << result << " J\n";
            break;
        default:
            std::cout << "Invalid Selection!\n";
    }
}
void newtonsLawofMotion() {
    int choice;
    double F, m, a, v, t, result;
    
    std::cout << "\n1. First Law of Newton (F = m * a)\n";
    std::cout << "2. Second Law of Newton (F = m * a)\n";
    std::cout << "3. Thirt Law of Newton (F1 = -F2)\n";
    std::cout << "4. Equation of motion.  (v = v0 + at)\n";
    std::cout << "Please Enter your Selection: ";
    std::cin >> choice;

    switch (choice) {
        // THE OTHER FORCE ACTING ON THE OBJECT AND KEEPING IT STATIONARY IS NEWTON'S 1ST LAW
        // THAT IS, IF A BOX OF "10KG" REMAINS STATIONARY ON THE GROUND, THE FORCE EXERTED BY THE GROUND ON THE BOX AND THE FORCE OF GRAVITY FROM ABOVE ARE EQUAL.
        case 1:
            std::cout << "DON'T FORGET!!! THE FIRST RULE AND THE SECOND RULE ARE ALMOST THE SAME..."<< std::endl;
            std::cout << "IF THE OBJECT IS STATIONARY OR MOVING AT A CONSTANT VELOCITY, THE FORCES ACTING ON THE OBJECT ARE EQUAL:"<< std::endl;
        
            std::cout << "Mass (kg): "; std::cin >> m;
            std::cout << "Angular Accelartion (m/s²): "; std::cin >> a;
            result = m * a;
            std::cout << "Newton: " << result << " N\n";
            break;
        case 2:
            std::cout << "Second Rule of Newtons Law:"<< std::endl;
            std::cout << "Mass (kg): "; std::cin >> m;
            std::cout << "Accelartion (m/s²): "; std::cin >> a;
            result = m * a;
            std::cout << "Newton: " << result << " N\n";
            break;
        case 3:
            std::cout << "First Force: "; std::cin >> F;
            std::cout << "Secon Force: "; std::cin >> F;
            result = -F;
            std::cout << "Newton: " << result << " N\n";
            break;
        case 4:
            std::cout << "First Speed (m/s): "; std::cin >> v;
            std::cout << "Time (speed): "; std::cin >> t;
            std::cout << "Angular Accelartion (m/s²): "; std::cin >> a;
            result = v + a * t;
            std::cout << "Speed: " << result << " m/s\n";
            break;
        default:
            std::cout << "Invalid Selection!\n";
    }
}
// ***KE MEAN KINETIC ENERJY***

void workEnergyandPower() {
    int choice;
    double F, d, t, W, P, E, m, h, g, result;
    std::cout << "\n1. Work (W = F * d)\n";
    std::cout << "2. Power (P = W / t)\n";
    std::cout << "3. Enerjy (E = F * d)\n";
    std::cout << "4. Potential Energy. (PE = m * g * h)\n";
    std::cout << "Please Enter your Selection: ";
    std::cin >> choice;

    switch (choice) {

        // F * d = work
    
        case 1:
            std::cout << "IF AN ANGLE IS GIVEN IN THE QUESTION, MULTIPLY THE RESULT BY COS. "<< std::endl;
            std::cout << "Newton (N): "; std::cin >> F;
            std::cout << "Distance (m): "; std::cin >> d;
            result = F * d;
            std::cout << "Power: " << result << " J\n";
            break;
        // W/t = power 
        case 2:
            std::cout << "Work (J): "; std::cin >> W;
            std::cout << "Time (second): "; std::cin >> t;
            result = W / t;
            std::cout << "Power: " << result << " W\n";
            break;
        // F * d looks same but totally different formuls.
        case 3:
            std::cout << "Newton (N): "; std::cin >> F;
            std::cout << "Distance (m): "; std::cin >> d;
            result = F * d;
            std::cout << "Enerjy: " << result << " J\n";
            break;
        // m * g * h gravtivical potential ...
        case 4:
            std::cout << "Mass (kg): "; std::cin >> m;
            std::cout << "Gravitical Accelartion (m/s²): "; std::cin >> g;
            std::cout << "Height (m): "; std::cin >> h;
            result = m * g * h;
            std::cout << "Potential Enerjy: " << result << " J\n";
            break;
        default:
            std::cout << "Invalid Selection!\n";
    }
}
void momentumandCollision() {
    int choice;
    double m1, m2, v1, v2, v1_prime, v2_prime, result;
    std::cout << "\n KE (KINETIC ENERJY)\n";
    std::cout << "\n1. Conservation of momentum.  (m1v1 + m2v2 = m1v1' + m2v2')\n";
    std::cout << "\n2.Elastic collision (m1v1 + m2v2 = m1v1' + m2v2' and CE is conserved). \n";
    std::cout << "\n3. Inelastic collision (there is CE loss or reduction).\n";
    std::cout << "Please Enter Your Selection: ";
    std::cin >> choice;

    switch (choice) {

        // Conservation of momentum in collision
        // Difference between total momentum before collision and total momentum after collision.
        
        case 1:
            std::cout << "First Stuff Mass (kg): "; std::cin >> m1;
            std::cout << "Second Stuff Mass (kg): "; std::cin >> m2;
            std::cout << "First Stuff Speed (m/s): "; std::cin >> v1;
            std::cout << "Second Stuff Speed (m/s): "; std::cin >> v2;
            std::cout << "New value of the velocity of the first object (m/s): "; std::cin >> v1_prime;
            std::cout << "New value of the velocity of the second object. (m/s): "; std::cin >> v2_prime;
            result = (m1 * v1 + m2 * v2) - (m1 * v1_prime + m2 * v2_prime);
            std::cout << "Momentum change is : " << result << " kg.m/s\n";
            break;
        case 2:
        // Change in kinetic energy (LOSS and CONSERVATION).
            std::cout << "First Stuff Mass (kg): "; std::cin >> m1;
            std::cout << "Second Stuff Mass (kg): "; std::cin >> m2;
            std::cout << "First Stuff Speed (m/s): "; std::cin >> v1;
            std::cout << "Second Stuff Speed (m/s): "; std::cin >> v2;
            std::cout << "New value of the velocity of the first object.(m/s): "; std::cin >> v1_prime;
            std::cout << "New value of the velocity of the second object. (m/s): "; std::cin >> v2_prime;
            result = 0.5 * m1 * (v1 * v1 - v1_prime * v1_prime) + 0.5 * m2 * (v2 * v2 - v2_prime * v2_prime);
            std::cout << "Kinetik enerji kaybı: " << result << " J\n";
            break;
        // Total energy change formula.
        case 3:
            std::cout << "First Stuff Mass (kg): "; std::cin >> m1;
            std::cout << "Second Stuff Mass (kg): "; std::cin >> m2;
            std::cout << "First Stuff Speed (m/s): "; std::cin >> v1;
            std::cout << "Second Stuff Speed (m/s): "; std::cin >> v2;
            std::cout << "New value of the velocity of the first object. (m/s): "; std::cin >> v1_prime;
            std::cout << "New value of the velocity of the second object. (m/s): "; std::cin >> v2_prime;
            result = 0.5 * m1 * (v1 * v1 + v2 * v2) - 0.5 * m1 * (v1_prime * v1_prime + v2_prime * v2_prime);
            std::cout << "Kinetic energy loss (inelastic collision).: " << result << " J\n";
            break;
        default:
            std::cout << "Invalid Selection!\n";
    }
}
void rotationalMotion (){
    int choice;
    double m, r, F, theta, v, I, omega, alpha, t, L, Ek, result;

    std::cout << "1. Angular Velocity Calculation. (angle/t)\n";
    std::cout << "2. Angular Acccelartion Calculation (rad/s - rad/s / t)\n";
    std::cout << "3. Torque Calculation (theta = theta * (pi / 180);)\n";
    std::cout << "4. Moment of Inertia Calculation.  (I = Kg * m^2(I))\n";
    std::cout << "5. Rotational Kinetic Energy Calculation.(1/2 * m * v^2 (KE))\n";
    std::cout << "6. Angular Moment Calculation (L)\n";
    std::cout << "Please Enter Your Selection: ";
    std::cin >> choice;

    switch (choice) {
        case 1:
            // Angular Speed Calculation
            std::cout << "Angular Displacement (rad): ";
            std::cin >> theta;
            std::cout << "Period (second): ";
            std::cin >> t;
            omega = theta / t;
            std::cout << "Angular Speed (omega): " << omega << " rad/s\n";
            break;

        case 2:
            // Angular accelartion calculation
            std::cout << "First Angular Speed (omega_i) (rad/s): ";
            std::cin >> omega;
            std::cout << "Last Angular Speed (omega_f) (rad/s): ";
            std::cin >> omega;
            std::cout << "Period (second): ";
            std::cin >> t;
            alpha = (omega - omega) / t;
            std::cout << "Angular Accelartion (alpha): " << alpha << " rad/s^2\n";
            break;

        case 3:
            // Torque Calculation
            std::cout << "Newton (N): ";
            std::cin >> F;
            std::cout << "Distance at Which the Force is Applied.(r) (metre): ";
            std::cin >> r;
            std::cout << "Angle Made by the Force with the lever (derece): ";
            std::cin >> theta;
            // convert angle to radian 
            theta = theta * (PI / 180);
            result = r * F * sin(theta);
            std::cout << "Tork (torque): " << result << " N.m\n";
            break;

        case 4:
            // Moment of Inertia
            std::cout << "Mass (kg): ";
            std::cin >> m;
            std::cout << "Radius (m): ";
            std::cin >> r;
            I = 0.5 * m * r * r;  // for the equal circle 
            std::cout << "Moment of Inertia (I): " << I << " kg.m^2\n";
            break;

        case 5:
            // Rotational kinetic energy calculation...
            std::cout << "Moment of inertia (I) (kg.m^2): ";
            std::cin >> I;
            std::cout << "Angular Velocity (omega) (rad/s): ";
            std::cin >> omega;
            Ek = 0.5 * I * omega * omega;
            std::cout << "Rotational kinetic energy. (Ek): " << Ek << " J\n";
            break;

        case 6:
            // Angular Momentum Calculation...
            std::cout <<"Moment of Inertia. (I) (kg.m^2): ";
            std::cin >> I;
            std::cout << "Angular Velocity (omega) (rad/s): ";
            std::cin >> omega;
            L = I * omega;
            std::cout << "Angular Momentum (L): " << L << " kg.m^2/s\n";
            break;

        default:
            std::cout << "Invalid Selection!\n";
            break;
    }
}
void angularMomentum(){

    int choice;
    double L, r, p, I, w, t, dL, dt, Loteleme, Ldonme, m, result;
    while (true) {
        std::cout << "Angular Momentum\n";
        std::cout << "1. Angular Momentum of the Particle (L = r * p * sin(theta))\n";
        std::cout << "2. Angular Momentum in Rotational Motion (L = I * w)\n";
        std::cout << "3. Torque ((tau) = dL / dt)\n";
        std::cout << "4. Angular Momentum Conservation (I1 * w1 = I2 * w2)\n";
        std::cout << "5. Postponence and Rotational Angular Momentum (L = L_Postponence + L_Rotational)\n";
        std::cout << "6. Exit\n";
        std::cout << "Please Enter Your Selection: ";
        std::cin >> choice;

        switch (choice) {
            case 1: {
                double theta;
                std::cout << "Radius (m): "; std::cin >> r;
                std::cout << "Mass (kg): "; std::cin >> m;
                std::cout << "Velocity (m/s): "; std::cin >> p;
                std::cout << "Theta (degree): "; std::cin >> theta;
                theta = theta * (PI / 180); 
                L = r * m * p * sin(theta);
                std::cout << "Angular Momentum: " << L << " kg·m²/s\n";
                break;
            }
            case 2: {
                std::cout << "Moment of Inertia (kg·m²): "; std::cin >> I;
                std::cout << "Angular Velocity (rad/s): "; std::cin >> w;
                L = I * w;
                std::cout << "Angular Momentum: " << L << " kg·m²/s\n";
                break;
            }
            case 3: {
                std::cout << "Differential Length: "; std::cin >> dL;
                std::cout << "Differential Time (s): "; std::cin >> dt;
                result = dL / dt;
                std::cout << "Torque (tau): " << result << " N·m\n";
                break;
            }
            case 4: {
                double I1, w1, I2, w2;
                std::cout << "Moment of Inertia (kg·m²): "; std::cin >> I1;
                std::cout << "Angular Velocity (rad/s): "; std::cin >> w1;
                std::cout << "Second Moment of Inertia (kg·m²): "; std::cin >> I2;
                w2 = (I1 * w1) / I2;
                std::cout << "New Angular Speed (w2): " << w2 << " rad/s\n";
                break;
            }
            case 5: {
                std::cout << "Postponence: "; std::cin >> Loteleme;
                std::cout << "Rotation: "; std::cin >> Ldonme;
                L = Loteleme + Ldonme;
                std::cout << "Total Angular Momentum: " << L << " kg·m²/s\n";
                break;
            }
            case 6:
                std::cout << "Logging Out...\n";
                return;
            default:
                std::cout << "Invalid Selection! Please Try Again...\n";
                }
            }
        }
void oscillationsandWaveMotion() {
    int choice; 
    double T, f, w, k, lambda, v, A, m, k_spring, x, E, U, K;
    
    while (true) {
        std::cout << "\nOscillations And Wave Motion\n";
        std::cout << "1. Period (T = 1 / f)\n";
        std::cout << "2. Frequency (f = 1 / T)\n";
        std::cout << "3. Angular Frequency (w = 2 * PI * f)\n";
        std::cout << "4. Wave Number (k = 2 * PI / lambda)\n";
        std::cout << "5. Wave Speed (v = f * lambda)\n";
        std::cout << "6. Simple Harmonic Motion Energy (E = 1/2 * k * A^2)\n";
        std::cout << "7. Kinetic Energy (K = 1/2 * m * v^2)\n";
        std::cout << "8. Potantial Energy (U = 1/2 * k * x^2)\n";
        std::cout << "9. Exit\n";
        std::cout << "Please Enter Your Selection: ";
        std::cin >> choice;
        
        switch (choice) {
            case 1:
                std::cout << "Frequency: (Hz): "; std::cin >> f;
                T = 1 / f;
                std::cout << "Period: " << T << " s\n";
                break;
            case 2:
                std::cout << "Period (s): "; std::cin >> T;
                f = 1 / T;
                std::cout << "Frequency: " << f << " Hz\n";
                break;
            case 3:
                std::cout << "Frequency (Hz): "; std::cin >> f;
                w = 2 * PI * f;
                std::cout << "Angular Frequency: " << w << " rad/s\n";
                break;
            case 4:
                std::cout << "Wave Height (m): "; std::cin >> lambda;
                k = 2 * PI / lambda;
                std::cout << "Wave Number: " << k << " 1/m\n";
                break;
            case 5:
                std::cout << "Frequency (Hz): "; std::cin >> f;
                std::cout << "Wave Height (m): "; std::cin >> lambda;
                v = f * lambda;
                std::cout << "Wave Speed: " << v << " m/s\n";
                break;
            case 6:
                std::cout << "Spring Constant (N/m): "; std::cin >> k_spring;
                std::cout << "Amplitude (m): "; std::cin >> A;
                E = 0.5 * k_spring * A * A;
                std::cout << "Total Energy: " << E << " J\n";
                break;
            case 7:
                std::cout << "Mass (kg): "; std::cin >> m;
                std::cout << "Velocity (m/s): "; std::cin >> v;
                K = 0.5 * m * v * v;
                std::cout << "Kinetic Energy: " << K << " J\n";
                break;
            case 8:
                std::cout << "Spring Constant (N/m): "; std::cin >> k_spring;
                std::cout << "Displacement (m): "; std::cin >> x;
                U = 0.5 * k_spring * x * x;
                std::cout << "Potential Energy: " << U << " J\n";
                break;
            case 9:
                std::cout << "Logging out...\n";
                return;
            default:
                std::cout << "Invalid Selection! Please Try Again.\n";
        }
    }
}
void fluidMechanics() {
 int choice; 
 double P, F, A, rho, h, g, Q, v, A1, A2, v1, v2, P1, P2, result;

    while (true) {
        std::cout << "\nFluid Mechanics\n";
        std::cout << "1. Pressure (P = F / A)\n";
        std::cout << "2. Hydrostatic Pressure (P = rho * g * h)\n";
        std::cout << "3. Debi (Q = A * v)\n";
        std::cout << "4. Bernoulli’s Principle (P1 + 0.5 * rho * v1^2 + rho * g * h1 = P2 + 0.5 * rho * v2^2 + rho * g * h2)\n";
        std::cout << "5. Exit\n";
        std::cout << "Please Enter Your Selection: ";
        std::cin >> choice;

        switch (choice) {
            case 1:
                std::cout << "Newton (N): "; std::cin >> F;
                std::cout << "Surface Area (m²): "; std::cin >> A;
                P = F / A;
                std::cout << "Pressure: " << P << " Pa\n";
                break;
            case 2:
                std::cout << "Fluid Density (kg/m³): "; std::cin >> rho;
                std::cout << "Gravitional Acceleration (m/s²): "; std::cin >> g;
                std::cout << "Height (m): "; std::cin >> h;
                P = rho * g * h;
                std::cout << "Hydrostatic Pressure: " << P << " Pa\n";
            case 3:
                std::cout << "Cross-sectional Area (m²): "; std::cin >> A;
                std::cout << "Flow Velocity (m/s): "; std::cin >> v;
                Q = A * v;
                std::cout << "Debi: " << Q << " m³/s\n";
                break;
            case 4:
                std::cout << "First Point Pressure (Pa): "; std::cin >> P1;
                std::cout << "First Point Velocity (m/s): "; std::cin >> v1;
                std::cout << "First Point Height (m): "; std::cin >> h;
                std::cout << "Second Point Pressure (Pa): "; std::cin >> P2;
                std::cout << "Second Point Velocity (m/s): "; std::cin >> v2;
                std::cout << "Second Point Height (m): "; std::cin >> h;
                result = (P1 + 0.5 * rho * v1 * v1 + rho * g * h) - (P2 + 0.5 * rho * v2 * v2 + rho * g * h);
                std::cout << "Bernoulli’s Principle Result: " << result << " Pa\n";
                break;
            case 5:
                std::cout << "Logging out...\n";
                return;
            default:
                std::cout << "Invalid Selection! Please Try Again.\n";
        }
    }
}
int main() {
    int choice;
    do {
        std::cout << "\n1. Motion and Dynamics \n";
        std::cout << "2. Newtons Law and Motion \n"; 
        std::cout << "3. Work Energy and Power \n";
        std::cout << "4. Momentum and Collision \n";
        std::cout << "5. Rotational Motion \n";
        std::cout << "6. Angular Momentum \n";
        std::cout << "7. Oscillations and Wavemotion \n";
        std::cout << "8. Fluid Mechanic \n";
        std::cout << "Please Enter Your Selection (0 exit): ";
        std::cin >> choice;
        
        switch (choice) {
            case 1:
                motionandDynamics();
                break;
            case 2:
                newtonsLawofMotion();
                break;
            case 3:
                workEnergyandPower();
                break;
            case 4:
                momentumandCollision();
                break;
            case 5:
                rotationalMotion();
                break;
            case 6:
                angularMomentum();
                break;
            case 7:
                oscillationsandWaveMotion();
                break;
            case 8:
                fluidMechanics();
                break;
            default:
                if (choice != 0) {
                    std::cout << "Invalid Selection!\n";
                }
                break;
        }
    } while (choice != 0);

    return 0;
}