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

// Random engine is declared as a static so that it can be used by many
// functions.
static default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
        
  num_particles = 75;

  // Normal distribution for sensor noise.
  normal_distribution<double> dist_x(0, std[0]);
  normal_distribution<double> dist_y(0, std[1]);
  normal_distribution<double> dist_theta(0, std[2]);  

  // Initialize all particles.
  for( int i = 0; i < num_particles; i++ ) {

    Particle p;
         
    p.id = i;
    p.x = x;
    p.y = y;
    p.theta = theta;
    p.weight = 1.0;

    // Add the Sensor noise
    p.x += dist_x(gen);
    p.y += dist_y(gen);
    p.theta += dist_theta(gen);

    particles.push_back(p);

  } // end of for loop

  // Initialization is completed.
  is_initialized = true;

} // end of ParticleFilter::init()


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // Normal distribution for sensor noise.
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  double div_vel_yawr = velocity / yaw_rate;
  double prod_yawr_delt = yaw_rate * delta_t;

  for( int i = 0; i < num_particles; i++ ) {

    // Handle yaw_rate close to 0.
    if ( fabs(yaw_rate) < 0.00001 ) {
      particles[i].x += velocity * delta_t * cos( particles[i].theta );
      particles[i].y += velocity * delta_t * sin( particles[i].theta );
    }
    else {

      particles[i].x += div_vel_yawr * ( sin(particles[i].theta + prod_yawr_delt) - sin(particles[i].theta) );
      particles[i].y += div_vel_yawr * ( cos(particles[i].theta) - cos(particles[i].theta + prod_yawr_delt) );
      particles[i].theta += prod_yawr_delt;
    }

    // Add noise.
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);


  } // end of for loop
    
} // end of ParticleFilter::prediction()

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.

  // For each observation, go through all the predictions and determine
  // which predicted landmark is nearest to the observed landmark.
  for (int i = 0; i < observations.size(); i++) {

    LandmarkObs observ = observations[i];

    double min_distance = numeric_limits<double>::max();
    int map_id = -1;

    for (int j = 0; j < predicted.size(); j++) {
      
      LandmarkObs pred = predicted[j];

      double obs_pred_distance = dist(observ.x, observ.y, pred.x, pred.y);

      if (obs_pred_distance < min_distance) {

        min_distance = obs_pred_distance;
        map_id = pred.id;
 
      } // end of if    

    } // end of for loop - predicted

    observations[i].id = map_id;    
    

  } // end of for loop - observations 

} // end of ParticleFilter::dataAssociation()


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	
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

  // Loop through all particles.
  for (int i = 0; i < num_particles; i++) {
    
    // Particle coordinates.
    double part_x = particles[i].x;
    double part_y = particles[i].y;
    double part_theta = particles[i].theta;

    vector<LandmarkObs> predictions;

    // Loop through all the landmarks.
    for (int j=0; j < map_landmarks.landmark_list.size(); j++) {

      // Map coordinates and ID of the landmark.
      float landm_x = map_landmarks.landmark_list[j].x_f;
      float landm_y = map_landmarks.landmark_list[j].y_f;
      int landm_id = map_landmarks.landmark_list[j].id_i;


      // Check if the landmark is within the sensor range.
      if (fabs(landm_x - part_x) <= sensor_range  &&
          fabs(landm_y - part_y) <= sensor_range ) {

        predictions.push_back( LandmarkObs{ landm_id, landm_x, landm_y } );

      } // end of if - sensor_range  
      

      // List of observations transformed from vehical coordinates to map
      // coordiantes.
      vector<LandmarkObs> transformed_obs;
      for (int k = 0; k < observations.size(); k++ ) {

        double trans_x = part_x + cos(part_theta) * observations[k].x - sin(part_theta) * observations[k].y;
        double trans_y = part_y + sin(part_theta) * observations[k].x + cos(part_theta) * observations[k].y; 
        transformed_obs.push_back( LandmarkObs{ observations[k].id, trans_x, trans_y } ); 

      } // end of for - observations.size


      // Data association for teh predictions and transformed observations for
      // the current particle.
      dataAssociation(predictions, transformed_obs);

      // Re-initialize the weight.
      particles[i].weight = 1.0;

      
      for (int l = 0; l < transformed_obs.size(); l++) {

        // Observation and associated prediction coordinates
        double obs_x, obs_y, pred_x, pred_y;
        obs_x = transformed_obs[l].x;
        obs_y = transformed_obs[l].y;
        int associated_prediction = transformed_obs[l].id;

        // The x and y coordinates of teh prediction associated with the current
        // observation.
        for (int m = 0; m < predictions.size(); m++) {

          if (predictions[m].id == associated_prediction) {

            pred_x = predictions[m].x;
            pred_y = predictions[m].y;

          } // end of if

        } // end of for - predictions.size()

        // Calculate the weight for this observation with multivariate
        // Gaussian.
        double s_x = std_landmark[0];
        double s_y = std_landmark[1];
        double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(pred_x -
obs_x,2)/(2*pow(s_x,2)) + (pow(pred_y - obs_y,2)/(2*pow(s_y,2))) ) );

        particles[i].weight *= obs_w;  
      } // end of for - transformed_obs.size()      

    } // end of for - landmark_list

  } // end of for - num_particles  

} // end of  ParticleFilter::updateWeights()


void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<Particle> new_particles;

  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {

    weights.push_back(particles[i].weight);

  } // end of for - num_particles


  uniform_int_distribution<int> uni_int_dist(0, num_particles-1);
  auto index = uni_int_dist(gen);

  double max_weight = *max_element(weights.begin(), weights.end());

  uniform_real_distribution<double> uni_real_dist(0.0, max_weight);

  double beta = 0.0;

  for (int i = 0; i < num_particles; i++ ) {
 
    beta += uni_real_dist(gen) * 2.0;

    while (beta > weights[index]) {

      beta -= weights[index];
      index = (index + 1) % num_particles;

    } // end of while

    new_particles.push_back(particles[index]);

  } // end of for

  particles = new_particles;

} // end of ParticleFilter::resample()


Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

} // end of


string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
} // end of 


string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
} // end of


string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
} // end of
