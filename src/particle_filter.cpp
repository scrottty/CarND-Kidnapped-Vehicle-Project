/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  
  // Create the normal distributions to create the particles with
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  // Set the number of particles
  num_particles = 10;
  
  // Create a random number generator
  default_random_engine gen;
  
  // Create num of particles with their positions initiliased around the GPS location
  for (int i=0; i<num_particles; i++)
  {
    Particle p;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles.push_back(p);
  }
  
  is_initialized = true;
  cout << "particle value, x: " << particles[0].x << " y: " << particles[0].y << " theta: " << particles[0].theta << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  // TEST VARIABLES
//  double x = 102;
//  double y = 65;
//  double theta = (5*M_PI)/8;
//  velocity = 110;
//  yaw_rate = M_PI/8;
  
  // Random number generator
  default_random_engine gen;
  
  for (int i=0; i<particles.size(); i++)
  {
    Particle p = particles[i];
    
    // Calc the change in position
    double x_change, y_change, theta_change;
    
    if (yaw_rate != 0)
    {
      x_change = velocity/yaw_rate * (sin(p.theta + yaw_rate*delta_t)-sin(p.theta));
      y_change = velocity/yaw_rate * (cos(p.theta)-cos(p.theta + yaw_rate*delta_t));
      theta_change = yaw_rate*delta_t;
    }
    else
    {
      x_change = velocity * cos(p.theta);
      y_change = velocity * sin(p.theta);
      theta_change = 0;
    }
    
    
//    cout << "x: " << x + x_change << " y: " << y + y_change << " theta: " << theta + theta_change << endl;
    
    // Create a normal distribution about the change
    normal_distribution<double> dist_x(x_change, std_pos[0]);
    normal_distribution<double> dist_y(y_change, std_pos[1]);
    normal_distribution<double> dist_theta(theta_change, std_pos[2]);
    
    // Update the particle positions with the noisey motion - predict new position
    p.x += dist_x(gen);
    p.y += dist_y(gen);
    p.theta += dist_theta(gen);
    
    particles[i] = p;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs>& observations, Map map_landmarks) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  
  // Loop through each observation
  for (int i=0; i<observations.size(); i++)
  {
    double minDist = 9999999;
    double minDist_idx = 0;
    
    // Loop through the landmarks and calculate the euclidean distance
    for (int j=0; j<map_landmarks.landmark_list.size(); j++)
    {
      double dist = sqrt(pow(observations[i].x - map_landmarks.landmark_list[j].x_f,2) + pow(observations[i].y - map_landmarks.landmark_list[j].y_f,2));
      
      // Find the minimum distance
      if (dist < minDist)
      {
        minDist = dist;
        minDist_idx = j;
      }
    }
    
    // The minimum distance index is the index of the closest landmarks
    // Set the id to the idx
    observations[i].id = minDist_idx;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  
  for (int i=0; i<particles.size(); i++)
  {
    Particle p = particles[i];
    std::vector<LandmarkObs> observations_shifted = observations;
    
    // Shift all of the observations to be as particle observations int the map coordinates
    for (int j=0; j< observations_shifted.size(); j++)
    {
      // TEST VARIABLES
//      double x_part = 4;
//      double y_part = 5;
//      double x_obs = 2;
//      double y_obs = 2;
//      double theta = -M_PI/2;
//      
//      p.x = x_part;
//      p.y = y_part;
//      p.theta = theta;
//      observations[j].x = x_obs;
//      observation s[j].y = y_obs;
      
      observations_shifted[j].x = p.x + (cos(p.theta)*observations[j].x) - (sin(p.theta) * observations[j].y);
      observations_shifted[j].y = p.y + (sin(p.theta)*observations[j].x) + (cos(p.theta) * observations[j].y);
      
//      cout << "xmap : " << observations_shifted[j].x << " ymap: " << observations_shifted[j].y << endl;
    }
    
    // Find the closet landmark to each of the points
    dataAssociation(observations_shifted, map_landmarks);
    
    double weight = 1;
    for (int j=0; j< observations_shifted.size(); j++)
    {
      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];
      double x = observations_shifted[j].x;
      double y = observations_shifted[j].y;
      double mu_x = map_landmarks.landmark_list[observations_shifted[j].id].x_f;
      double mu_y = map_landmarks.landmark_list[observations_shifted[j].id].y_f;
      
      double normaliser = 1/(2*M_PI*sig_x*sig_y);
      
      double exponent = -(pow(x-mu_x,2)/(2*pow(sig_x,2)) + pow(y-mu_y,2)/(2*pow(sig_y,2)));
      
      weight *= normaliser * pow(M_E,exponent);
    }
    p.weight = weight;
    particles[i] = p;
  }
  
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  // New vector for the resampled particles
  std::vector<Particle> particles_resampled;
  
  // Make a vector with all of the weights
  std::vector<double> weights;
  for (int i=0; i<particles.size(); i++)
  {
    weights.push_back(particles[i].weight);
  }
  
  // distribution with to select from and random generator
  std::discrete_distribution<double> dist (weights.begin(), weights.end());
  std::default_random_engine gen;
  
  // Select new particles
  // Add particle to vector selected from weights in dist
  for (int i=0; i<particles.size(); i++)
  {
    particles_resampled.push_back(particles[dist(gen)]);
  }
  
  // Replace particles with the resampled one
  particles = particles_resampled;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
