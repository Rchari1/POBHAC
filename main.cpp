#include <iostream>
#include <cmath>



// Global Variables
const double G = 6.67430e-11;  // Gravitational constant
const double c = 2.998e8;      // Speed of light
const double epsilon = 0.1;   // Radiative efficiency for a Schwarzschild black hole
const double sigma_T = 6.652e-29; // Thomson scattering cross-section
const double sigma = 1e-26;   // Cross-section for the spallation reaction near a black hole
const double f_p = 1e16;      // Proton flux
const double R = 1e12;        // Radius of the accretion disk
const double Omega = 1e-7;    // Angular velocity of the accretion disk
double n_O = 1e14;            // Initial number density of oxygen
double n_Be = 1e13;           // Initial number density of beryllium
const double M_J = 1.898e27;  // Mass of Jupiter
const double T = 1000;        // Temperature
const double m = 1e-25;       // Mean mass of elements in dust particles
const double m_b = 16;     // Mass of beryllium
const double m_u = 9.012;  // Atomic mass unit
const double k = 1.380649e23; // Boltzman constant
double Mdot = 0.0;
double totalMassLost = 0.0;
double f_pl = 0.01; 
double M1 = 1.989e+38; // Mass of the Black Hole
double M2 = 5.972e24; // Mass of planet
double M3 = 7.342e22; // Mass of exomoon 

// Define threshold for "consumed" status
const double CONSUMED_THRESHOLD = 0.05;  // 5% of original mass

// Timestepping
const double dt_max = 1.0;
const double dt_min = 0.01;
const double critical_distance = 1e8;

struct Body {
    double mass;
    double radius;
    double x, y;
    double vx, vy;
};

void computeForce(const Body& a, const Body& b, double& fx, double& fy) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double r = std::sqrt(dx*dx + dy*dy);
    double F = G * a.mass * b.mass / (r*r*r);
    fx = F * dx;
    fy = F * dy;
}

bool isTidallyDisrupted(const Body& b, const Body& blackHole) {
    double r = std::sqrt((b.x - blackHole.x)*(b.x - blackHole.x) + (b.y - blackHole.y)*(b.y - blackHole.y));
    double tidalRadius = (b.radius) * std::pow((blackHole.mass / b.mass), 1.0/3.0);
    return r < tidalRadius;
}
/*
double calculate_L_dot(const Body& blackHole, const Body& disruptedBody) {
    double dx = disruptedBody.x - blackHole.x;
    double dy = disruptedBody.y - blackHole.y;
    double relative_velocity = std::sqrt(std::pow(disruptedBody.vx - blackHole.vx, 2) + std::pow(disruptedBody.vy - blackHole.vy, 2));
    return blackHole.mass * disruptedBody.mass * relative_velocity / (dx*dx + dy*dy);
}


void accretionFromAngularMomentum(Body& blackHole, const Body& disruptedBody) {
    double L_dot = calculate_L_dot(blackHole, disruptedBody);
    Mdot = L_dot / (2 * M_PI * R * R * Omega); 
    double accretedMass = Mdot * 1.0;
    blackHole.mass += std::min(accretedMass, disruptedBody.mass);
}

*/


double accretionFromAngularMomentum(Body& blackHole, const Body& disruptedBody) {
    double L_dot = ((4 * M_PI * G * M1 * c)/ sigma_T);
    Mdot = L_dot / (epsilon * c * c); 
    double accretedMass = Mdot * 1.0;
    double actualAccretedMass = std::min(accretedMass, disruptedBody.mass);
    blackHole.mass += actualAccretedMass;
    return actualAccretedMass;
}



double calculate_mixing_ratio(double f_pl, double Mdot,  Body& blackHole) {
    double blackHole_mass_J = blackHole.mass / M_J; // Black hole mass in Jupiter masses
    double P = 1e-6 * 1e5; // Pressure in Pascals (1 microbar)
    double n_b = P / (k * T); // Background number density in m^-3 (k is the Boltzmann constant)
    
    // Directly using Mdot from the function arguments
    double F = -6.5e9 * f_pl * (epsilon * Mdot / 1e8) * std::pow(R / 1737e3, -2); // Particle flux

    double n_i_over_n_b = 4.1e-3 * f_pl * (epsilon * Mdot / 1e8) * std::pow(blackHole_mass_J, -1) * std::sqrt(T / 1000) * std::sqrt((m_b * m_u) / (m * (m + m_b)));

    n_i_over_n_b /= n_b;

    return n_i_over_n_b;
}


int main() {
    Body blackHole = {M1, 1e12, 0, 0, 0, 0};
    Body planet = {M2, 6371e6, 1e9, 0, 0, 2e4};
    Body exomoon = {M3, 1737e3, 1e9 + 1e8, 0, 0, 2e4 + 1e3};

    double dt = 1.0;

    for (int step = 0; step < 10000; ++step) {

    double distanceExomoonToBlackHole = std::sqrt((exomoon.x - blackHole.x)*(exomoon.x - blackHole.x) + (exomoon.y - blackHole.y)*(exomoon.y - blackHole.y));

    // Variable time-stepping
    if (distanceExomoonToBlackHole < critical_distance) {
        dt = dt_min;
    } else {
        dt = dt_max;
    }

    if (isTidallyDisrupted(planet, blackHole)) {
        double accretedFromPlanet = accretionFromAngularMomentum(blackHole, planet);
        planet.mass -= accretedFromPlanet;
        totalMassLost += accretedFromPlanet;
    }

    if (isTidallyDisrupted(exomoon, blackHole)) {
        double accretedFromExomoon = accretionFromAngularMomentum(blackHole, exomoon);
        exomoon.mass -= accretedFromExomoon;
        totalMassLost += accretedFromExomoon;
        double delta_n_Be = sigma * f_p * n_O * dt;
        n_Be += delta_n_Be;
        n_O -= delta_n_Be;
    }




        double fpx, fpy, fmx, fmy;
        computeForce(planet, blackHole, fpx, fpy);
        computeForce(exomoon, blackHole, fmx, fmy);

        planet.vx += fpx / planet.mass * dt;
        planet.vy += fpy / planet.mass * dt;
        exomoon.vx += fmx / exomoon.mass * dt;
        exomoon.vy += fmy / exomoon.mass * dt;

        planet.x += planet.vx * dt;
        planet.y += planet.vy * dt;
        exomoon.x += exomoon.vx * dt;
        exomoon.y += exomoon.vy * dt;

        // Mixing ratio

        double mixing_ratio = calculate_mixing_ratio(f_pl, Mdot, blackHole);

        // Check if the planet or exomoon is consumed
        if (planet.mass < M2 * CONSUMED_THRESHOLD) {
            std::cout << "Planet effectively consumed by the black hole at step " << step << "." << std::endl;
            break;  // Exit the simulation loop
        }

        if (exomoon.mass < M3 * CONSUMED_THRESHOLD) {
            std::cout << "Exomoon effectively consumed by the black hole at step " << step << "." << std::endl;
            break;  // Exit the simulation loop
        }


         // Output
        if (step % 100 == 0) {
            double distanceToBlackHole = std::sqrt((exomoon.x - blackHole.x)*(exomoon.x - blackHole.x) + (exomoon.y - blackHole.y)*(exomoon.y - blackHole.y));



            std::cout << "Step: " << step
                    << " Planet: (" << planet.x << ", " << planet.y << ")"
                    << " Exomoon: (" << exomoon.x << ", " << exomoon.y << ")"
                    << " Beryllium to Oxygen ratio: " << n_Be / n_O 
                    << " Distance (Exomoon to Black Hole): " << distanceToBlackHole 
                    << " Mixing Ratio: " << mixing_ratio
                    << "Total Mass Lost: " << totalMassLost 
                    << std::endl;
}


    }

/*
Simulating the interaction between a planet, its exomoon, and a black hole. This includes the gravitational dynamics, 
tidal disruptions, and the spallation reactions leading to beryllium production as the exomoon gets close to the black hole.
*/

    return 0;
}
