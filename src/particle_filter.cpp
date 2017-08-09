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
#include<limits>


using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//Create random normal distribution

	default_random_engine gen;
	normal_distribution<double> dist(0.0,1.0);

	// Check if dim=3 or the length = 3 ??  Guessing length = 3
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	// Set this as startng value
	num_particles = 50;

	// Initialize all the particle weights to 1
	weights.resize(num_particles);
	fill(weights.begin(),weights.end(),1) ;

	// Set all the initial values for the paricles
	particles.resize(num_particles);

	for(unsigned int i =0; i < num_particles; i++){
		particles[i].id = i ;
		particles[i].x = x + std_x * dist(gen);
		particles[i].y = y + std_y * dist(gen);
		particles[i].theta = theta + std_theta * dist(gen);
		particles[i].weight = 1.0 ;
	}

	// set flag is_initialized
	is_initialized = true;


}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//Create random normal distribution
	default_random_engine gen;
	normal_distribution<double> dist(0.0,1.0);


	double std_x = std_pos[0] ;
	double std_y = std_pos[1] ;
	double std_theta = std_pos[2];

	// Update each of the particles using the bicycle model

	for(unsigned int i =0; i < num_particles; i++){
		double theta = particles[i].theta ;

		// non-zero yaw_rate
		if(fabs(yaw_rate) > 0.001){

			particles[i].x +=  (velocity / yaw_rate) * ( sin(theta + yaw_rate * delta_t) - sin(theta) ) +
					std_x * dist(gen);
			particles[i].y += (velocity / yaw_rate) * ( cos(theta) - cos(theta + yaw_rate * delta_t)) +
					std_y * dist(gen);
			particles[i].theta += yaw_rate *delta_t + std_theta * dist(gen);
		}
		// use simple model when yaw_rate is ~zero
		else {

			particles[i].x += velocity * cos(theta) * delta_t + std_x * dist(gen) ;
			particles[i].y += velocity * sin(theta) * delta_t + std_y * dist(gen) ;
			particles[i].theta += yaw_rate * delta_t + std_theta * dist(gen);

		}
	}

}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	double distance;
	double obs_size = observations.size() ;
	double pred_size = predicted.size() ;

	// Nearest Neighbour - Is there a better way to do it than O(mn) ?
		for (unsigned int i = 0; i < obs_size; i++){

			double minimum_dist = std::numeric_limits<double>::max() ;
			int pred_index = -1;

			for (int j=0;j < pred_size; j++){

				distance = dist(observations[i].x, observations[i].y,
											predicted[j].x, predicted[j].y);
				if (distance < minimum_dist){
					minimum_dist = distance;
					pred_index = j;   // Provide the index of the predicted measurement
				}
			}
			observations[i].id  = pred_index; // Keep track of Array index
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
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

    double std_x = std_landmark[0];

    double std_y = std_landmark[1];

    double sum_of_weights = 0.0 ;

    unsigned int num_obs = observations.size();

     weights.clear() ;

    //Iterate over every particle
    for(unsigned int i = 0; i < num_particles; i++){ // i variable to iterate over particle s

		// Create a vector of observations with the MAP co-ordinates as references
		vector<LandmarkObs> map_obs;
		map_obs.resize(num_obs) ;

		// Store the theta value that will be used repeatedly
		double theta = particles[i].theta ;

		// Rotation + Translation transformation for every observation
		for(unsigned int j = 0; j < num_obs; j ++ ){ // j variable to iterate over observations

			map_obs[j].id = -1; // This has to be matched in the dataAssociation Step

			map_obs[j].x = observations[j].x * cos(theta) -( observations[j].y * sin(theta) ) + particles[i].x ;

			map_obs[j].y = observations[j].x * sin(theta) + observations[j].y * cos(theta) + particles[i].y ;

		}

		vector<LandmarkObs> predicted ;

		for ( unsigned int j = 0; j < map_landmarks.landmark_list.size(); j ++){
					double distance = dist(particles[j].x, particles[j].y,
										   map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) ;

			// Check if in sensor range
			if(distance < sensor_range){
				LandmarkObs temp;
				temp.id = map_landmarks.landmark_list[j].id_i ;
				temp.x = map_landmarks.landmark_list[j].x_f ;
				temp.y = map_landmarks.landmark_list[j].y_f ;

				predicted.push_back(temp) ;
			}


		}

		// Match the observations with the landmarks based on predictions
		dataAssociation(predicted, map_obs);

      //Calculate Bivariate Probability Product and Update the weight values
		double prob = 1 ;

		for(unsigned int j = 0; j < num_obs; j++){
			int predicted_id =  map_obs[j].id ; // This is the array position of the predicted id
			if(predicted_id > -1){ // match was found
				prob *= bivariateGaussian(predicted[predicted_id].x, map_obs[j].x,
										  predicted[predicted_id].y, map_obs[j].y, std_x, std_y);
			}
		}

		weights.push_back(prob);
		particles[i].weight = prob;
		sum_of_weights += prob ;
	}

    // Normalize all the weights

    for(unsigned int k =0; k < num_particles; k++){
    	particles[k].weight /= sum_of_weights ;
    	weights[k] = particles[k].weight ;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;

	discrete_distribution<> dist_particles(weights.begin(), weights.end());

    vector<Particle> resampled_particles;

    resampled_particles.resize(num_particles);

    // Resample
    for(int i = 0; i < num_particles; i++){

    	resampled_particles[i] = particles[dist_particles(gen)];

    }
    particles = resampled_particles;

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
