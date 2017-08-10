#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include <random>
#include <cstdio>
#include <cstdlib>
#include "spline.h"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
// DEBUG
ofstream myfile;

// constants

// saved lists
vector<double> car_accel_list;
vector<double> car_velo_list;
vector<double> car_pos_list;
vector<double> car_dpos_list;
string current_state = "KL";
// booleans
double START_UP = 1;
int lane = 1;
bool warning = false;
double SPACING = 30;
// constants
//double SPEED_LIMIT = 18.5;  // previously 19 can be partially working  3.6 miles best
double SPEED_LIMIT = 21;  // previously 19 can be partially working  3.6 miles best
double BUFFER_ZONE = 30.0;
double BUFFER_ZONE2 = 60.0;
double FREE_ZONE = 100.0;
double BUBBLE_DIST = 20.0; // change from 10 to 20
double max_s = 6945.554;
double v_targ = 0;

// weiths
double SPEED_WEIGHT = 1.0;
double DISTANCE_WEIGHT = 2.0;
double LANE_WEIGHT = 0.5;
double BUBBLE_DIST_WEIGHT = 99.0; // change from 10 to 20

// for convenience
using json = nlohmann::json;

// state flags
double TURN = 1;
double LTURN = 1;
double RTURN = 1;
double NOT_TURN = 0;
double TURN_COOLDOWN = 5;
double KEEP_LANE = 1;
double GOING_LEFT = 0;
double GOING_RIGHT = 0;


// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2) {
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y) {
	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++) {
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen) {
			closestLen = dist;
			closestWaypoint = i;
		}
	}
	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4) {
		closestWaypoint++;
	}
  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0) {
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point
	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef) {
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++) {
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y) {
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) )) {
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};
}

// logitic function
double logistic(double x) {
  return (2.0 / (1 + exp(-x)) - 1.0);
}

/*
double wrap_distance(double s_input) {
    return s_input - max_s*floor(s_input/max_s);
} */

// calculate a polynomial equation
double calc_equation(vector<double> coefficients, double T) {
  /*
  Takes the coefficients of a polynomial and creates a function of
  time from them.
  */
  double total = 0.0;
  for (double i=0; i<coefficients.size(); i++) {
    total += coefficients[i]*pow(T,i);
  }
  return total;
}

vector<double> differentiate(vector<double> coefficients){
  /*
  Calculates the derivative of a polynomial and returns
  the corresponding coefficients.
  */
  vector<double> new_cos;
  new_cos.clear();
  for (double i = 1.0; i<coefficients.size(); i++) {
    new_cos.push_back(i*coefficients[i]);
  }
  return new_cos;
}

// Begin the COST Function designs
// a version from python udacity implementation


int check_lane (double d_val) {
  if (d_val > 0.0 && d_val <= 4.0) {
    return 0; // Lane 0 is left most lane
  } else if (d_val > 4.0 && d_val <= 8.0) {
    return 1; // Lane 1 is middle lane
  } else if (d_val > 8.0 && d_val < 12.0) {
    return 2; // Lane 2 is right most lane
  } else {
    return -1;
  }
}

// lane middle definitions
double find_lane_mid (int lane) {
  if (lane == 0) {
    return 2.0; // Lane 0 is left most lane
  } else if (lane == 1) {
    return 6.0; // Lane 1 is middle lane
  } else if (lane == 2) {
    return 10.0; // Lane 2 is right most lane there's some skew near the end
  } else {
    return 0;
  }
}

vector<double> JMT(vector<double> start, vector <double> end, double T) {
    MatrixXd A = MatrixXd(3, 3);
    double T_2 = T*T;
    double T_3 = T_2*T;
    double T_4 = T_2*T_2;
    double T_5 = T_3*T_2;

	A <<   T_3,   T_4,    T_5,
			 3*T_2, 4*T_3,  5*T_4,
			 6*T,  12*T_2, 20*T_3;

	MatrixXd B = MatrixXd(3,1);
	B << end[0]-(start[0] + start[1] * T + .5*start[2]*T_2),
			 end[1]-(start[1] + start[2] * T),
			 end[2]- start[2];

	MatrixXd Ai = A.inverse();

	MatrixXd C = Ai*B;

	vector <double> result = {start[0], start[1], 0.5*start[2]};
	for(int i = 0; i < C.size(); i++) {
	    result.push_back(C.data()[i]);
	}
  return result;
}

// AT cars
//cout << " -------------------------------- " << endl;
//cout << "id of car: " << sensor_fusion[i][0] << endl;
//cout << "x of car: " << sensor_fusion[i][1] << endl;
//cout << "y of car: " << sensor_fusion[i][2] << endl;
//cout << "vx of car: " << sensor_fusion[i][3] << endl;
//cout << "vy of car: " << sensor_fusion[i][4] << endl;
//cout << "s of car: " << sensor_fusion[i][5] << endl;
//cout << "d of car: " << sensor_fusion[i][6] << endl;
double calc_speed(vector<double> other_car) {
  double speed = sqrt(other_car[3]*other_car[3]+other_car[4]*other_car[4]);
  return speed;
}

double s_dist_away(vector<double> other_car, double my_car_s) {
  double distance_away = other_car[5] - my_car_s;
  return distance_away;
}

double calc_dist(double other_car_s, double my_car_s) {
  double distance_away = other_car_s - my_car_s;
  return distance_away;
}

vector<double> check_car_front(vector<vector<double>> sensors, int my_car_lane, double my_car_s) {
  vector<vector<double>> cars_in_same_lane;
  vector<double> target_car = {-1.0};
  cars_in_same_lane.clear();

  for (int i=0; i<sensors.size(); i++) {
    int targ_car_lane = check_lane(sensors[i][6]);
    if (targ_car_lane == my_car_lane) {
      cars_in_same_lane.push_back(sensors[i]);
    }
  }
  cout << "# of cars in same lane: " << cars_in_same_lane.size() << endl;
  double max_dist = 999999.0;
  for (int cr=0; cr < cars_in_same_lane.size(); cr++) {
    //cout << "here are the cars in my lane: " << cars_in_same_lane[cr][0] << endl;
    //cout << "here are the cars s dist: " << cars_in_same_lane[cr][5] << endl;
    // wrapping the numbers
    double curr_dist = cars_in_same_lane[cr][5] - my_car_s;
    if (curr_dist < 0) {
      curr_dist += max_s;
    }

    if (curr_dist > 0 && curr_dist < max_dist) {
      target_car = cars_in_same_lane[cr];
      max_dist = curr_dist;
    }
  }
  cout << "car s: " << target_car[0] << endl;
  return target_car;
}

vector<double> check_car_left_front(vector<vector<double>> sensors, int my_car_lane, double my_car_s) {
  vector<vector<double>> cars_in_left_lane;
  vector<double> target_car = {-1.0};
  cars_in_left_lane.clear();

  for (int i=0; i<sensors.size(); i++) {
    int targ_car_lane = check_lane(sensors[i][6]);
    if (targ_car_lane == (my_car_lane-1)) {
      cars_in_left_lane.push_back(sensors[i]);
    }
  }

  double max_dist = 99999;

  for (int cr=0; cr < cars_in_left_lane.size(); cr++) {
    // wrapping the numbers
    double curr_dist = cars_in_left_lane[cr][5] - my_car_s;
    if (curr_dist < 0) {
      curr_dist += max_s;
    }
    if (curr_dist < max_dist && curr_dist >0) {
      target_car = cars_in_left_lane[cr];
      max_dist = curr_dist;
    }
  }
  return target_car;
}

vector<double> check_car_left_rear(vector<vector<double>> sensors, int my_car_lane, double my_car_s) {
  vector<vector<double>> cars_in_left_lane;
  vector<double> target_car = {-1.0};
  cars_in_left_lane.clear();

  for (int i=0; i<sensors.size(); i++) {
    int targ_car_lane = check_lane(sensors[i][6]);
    if (targ_car_lane == (my_car_lane-1)) {
      cars_in_left_lane.push_back(sensors[i]);
    }
  }

  double max_dist = -99999;

  for (int cr=0; cr < cars_in_left_lane.size(); cr++) {
    double curr_dist = my_car_s - cars_in_left_lane[cr][5];
    if (curr_dist < 0) {
      curr_dist += max_s;
      curr_dist *= -1;
    }
    if ((curr_dist > max_dist) && (curr_dist < 0)) {
      target_car = cars_in_left_lane[cr];
      max_dist = curr_dist;
    }
  }
  return target_car;
}

vector<double> check_car_right_front(vector<vector<double>> sensors, int my_car_lane, double my_car_s) {
  vector<vector<double>> cars_in_right_lane;
  vector<double> target_car = {-1.0};
  cars_in_right_lane.clear();

  for (int i=0; i<sensors.size(); i++) {
    int targ_car_lane = check_lane(sensors[i][6]);
    if (targ_car_lane == (my_car_lane+1)) {
      cars_in_right_lane.push_back(sensors[i]);
    }
  }

  double max_dist = 99999;

  for (int cr=0; cr < cars_in_right_lane.size(); cr++) {
    double curr_dist = cars_in_right_lane[cr][5] - my_car_s;
    if (curr_dist < 0) {
      curr_dist += max_s;
    }
    if (curr_dist < max_dist && curr_dist >0) {
      target_car = cars_in_right_lane[cr];
      max_dist = curr_dist;
    }
  }
  return target_car;
}

vector<double> check_car_right_rear(vector<vector<double>> sensors, int my_car_lane, double my_car_s) {
  vector<vector<double>> cars_in_right_lane;
  vector<double> target_car = {-1.0};
  cars_in_right_lane.clear();

  for (int i=0; i<sensors.size(); i++) {
    int targ_car_lane = check_lane(sensors[i][6]);
    if (targ_car_lane == (my_car_lane+1)) {
      cars_in_right_lane.push_back(sensors[i]);
    }
  }

  double max_dist = -99999;

  for (int cr=0; cr < cars_in_right_lane.size(); cr++) {
    double curr_dist = my_car_s - cars_in_right_lane[cr][5];
    if (curr_dist < 0) {
      curr_dist += max_s;
      curr_dist *= -1;
    }
    if ((curr_dist > max_dist) && (curr_dist < 0)) {
      target_car = cars_in_right_lane[cr];
      max_dist = curr_dist;
    }
  }
  return target_car;
}


double calc_cost(double my_car_s, double my_car_speed, int my_lane, vector<double> front_car, vector<double> back_car, double left_move, double right_move, vector<vector<double>> predictions, int prev_path_size, double time_inc) {
  double total_cost = 0;
  // speed cost
  double front_car_speed;
  double back_car_speed;
  double back_car_dist;
  double front_car_dist;

  // find out the max_dist to normalize for distance cost
  /*
  double max_dist = 0;

  for (auto v : predictions) {
    double d_away = abs(s_dist_away(v, my_car_s));
    if (d_away > max_dist) {
      max_dist = d_away ;
    }
  } */


  if (front_car[0] != -1) {
    double front_car_s = front_car[5];
    double front_car_speed = calc_speed(front_car);
    // predict future car position
    double front_car_future_s = front_car_s + (double)prev_path_size*time_inc*front_car_speed;
    double front_car_dist = calc_dist(front_car_future_s, my_car_s);


    /*
    if (front_car_speed < SPEED_LIMIT) {
      total_cost += 1;
    }
    if (front_car_speed > my_car_speed) {
      total_cost += 1/(1+abs(front_car_speed-my_car_speed));
    } else {
      total_cost += 1;
    } */
    /*
    if (front_car_dist <= 10.0) {
      total_cost += 5;
    } else if (front_car_dist > 10.0 && front_car_dist<= 20.0) {
      total_cost += 4;
    } else if (front_car_dist > 20.0 && front_car_dist<= 30.0) {
      total_cost += 3;
    } else if (front_car_dist > 30.0 && front_car_dist<= 40.0) {
      total_cost += 2;
    } else if (front_car_dist > 40.0 && front_car_dist<= 50.0) {
      total_cost += 1;
    } else {
      total_cost += 0;
    }
  } */
    // lane costs
    if ((left_move ==1) || (right_move==1)) {
      total_cost += LANE_WEIGHT;
    }
    //cout << "total cost after lane weight: " << total_cost << endl;
    if (front_car_speed > SPEED_LIMIT) {
      total_cost += 0;
    } else {
      total_cost += SPEED_WEIGHT*((SPEED_LIMIT-front_car_speed)/SPEED_LIMIT);
    }
    //cout << "total cost after speed weight: " << total_cost << endl;
    double FAR_SPACE = 120.0;
    if (front_car_dist > FAR_SPACE) {
      total_cost += 0;
    } else if (front_car_dist <= FAR_SPACE && front_car_dist > 0){
      total_cost += DISTANCE_WEIGHT*((FAR_SPACE-front_car_dist)/FAR_SPACE);
    } else {
      total_cost += 0;
    }
    //cout << "total cost after after_dist weight: " << total_cost << endl;
  }
  /*
  if (back_car[0] != -1) {
    back_car_speed = calc_speed(back_car);
    back_car_dist = s_dist_away(back_car, my_car);

    if (abs(back_car_dist) <= 10.0) {
      total_cost += 5;
    } else if (abs(back_car_dist) > 10.0 && abs(back_car_dist)<= 20.0) {
      total_cost += 4;
    } else if (abs(back_car_dist) > 20.0 && abs(back_car_dist)<= 30.0) {
      total_cost += 3;
    } else if (abs(back_car_dist) > 30.0 && abs(back_car_dist)<= 40.0) {
      total_cost += 2;
    } else if (abs(back_car_dist) > 40.0 && abs(back_car_dist)<= 50.0) {
      total_cost += 1;
    } else {
      total_cost += 0;
    }
  } */


  //total_cost += SPEED_WEIGHT*(SPEED_LIMIT-my_car[1])*(SPEED_LIMIT-my_car[1]);
  // distance cost

  // closer the distance the higher the cost
  // farther the distance the lower the cost
  //total_cost += DISTANCE_WEIGHT/(1+abs(front_car_dist-my_car));

  // penalize multiple lane change

  if (left_move == 1 || right_move == 1) {
    for (auto v : predictions) {
      double v_s = v[5];
      double v_lane = check_lane(v[6]);
      double v_speed = calc_speed(v);
      // predict future car position
      double v_future_s = v_s + (double)prev_path_size*time_inc*v_speed;
      double v_dist = abs(calc_dist(v_future_s, my_car_s));

      if ( (v_lane+1 == my_lane) && left_move == 1 && v_dist < BUBBLE_DIST) {
        //cout << "I shoudln't turn left" << endl;
        total_cost += 99.0;
        //cout << "this is total cost in not turn left: " << total_cost << endl;

      } else if ( (v_lane-1 == my_lane) && right_move == 1 && v_dist < BUBBLE_DIST) {
        //cout << "I shoudln't turn right" << endl;
        total_cost += 99.0;
        //cout << "this is total cost in not turn right: " << total_cost << endl;

      } else {
        total_cost += 0;
      }
    }
  }
  //cout << "this is total cost inside: " << total_cost << endl;

  return total_cost;
}

void update_velocity(double car_s, vector<double> front_car, int prev_path_size, double time_inc) {
  if (front_car[0] != -1) {
    double front_car_s = front_car[5];
    double front_car_speed = calc_speed(front_car);
    // predict future car position
    double front_car_future_s = front_car_s + (double)prev_path_size*time_inc*front_car_speed;
    double front_car_dist = calc_dist(front_car_future_s, car_s);

    cout << "here is front car speed: " << front_car_speed << endl;
    cout << "here is front car future s: " << front_car_future_s << endl;
    if (front_car_dist >0 && (front_car_dist)<BUFFER_ZONE) {
      warning = true;
    } else {
      warning = false;
    }
    cout << "here is my target velocity: " << v_targ << endl;
  } else {
    warning = false;
  }

  if (warning) {
    v_targ -= 0.1;
  } else if (v_targ < SPEED_LIMIT-.1) {
    v_targ += 0.1;
  }
}


int main() {
  uWS::Hub h;

  // initialization
  int first_start = 1;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }


  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
            int car_lane = check_lane(car_d);

            auto delta = {-10.0, 0.0, 0.0, 0.0, 0.0 ,0.0};

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

            // attempting walkthrough version

            // we will create some evenly spaced x,y waypoints at 30m
            // then interpolate with spline
            vector <double> x_val;
            vector <double> y_val;
            double time_inc = 0.02;

            double start_x = car_x;
            double start_y = car_y;
            double start_yaw = deg2rad(car_yaw);

            int prev_path_size = previous_path_x.size();

            if (prev_path_size>0) {
              car_s = end_path_s;
            }

            vector<double> front_car = check_car_front(sensor_fusion, car_lane, car_s);
            vector<double> back_car = {-1};
            vector<double> left_front_car = check_car_left_front(sensor_fusion, car_lane, car_s);
            vector<double> left_back_car = check_car_left_rear(sensor_fusion, car_lane, car_s);
            vector<double> right_front_car = check_car_right_front(sensor_fusion, car_lane, car_s);
            vector<double> right_back_car = check_car_right_rear(sensor_fusion, car_lane, car_s);

            double front_car_speed = calc_speed(front_car);
            double left_front_car_speed = calc_speed(left_front_car);
            double right_front_car_speed = calc_speed(right_front_car);

            double front_car_dist = s_dist_away(front_car, car_s);
            double left_front_car_dist = s_dist_away(left_front_car, car_s);
            double right_front_car_dist = s_dist_away(right_front_car, car_s);

            cout << "front car dist: " <<  front_car_dist << endl;
            cout << "L front car dist: " <<  left_front_car_dist << endl;
            cout << "R front car dist: " <<  right_front_car_dist << endl;
            cout << "front car speed: " <<  front_car_speed << endl;
            cout << "L front car speed: " <<  left_front_car_speed << endl;
            cout << "R front car speed: " <<  right_front_car_speed << endl;

            double front_cost = 99999;
            double left_cost = 99999;
            double right_cost = 99999;

            if (car_lane ==1) {
              front_cost = calc_cost(car_s, car_speed*0.44704, car_lane, front_car, back_car, NOT_TURN, NOT_TURN, sensor_fusion, prev_path_size, time_inc);
              right_cost = calc_cost(car_s, car_speed*0.44704, car_lane, right_front_car, right_back_car, NOT_TURN, RTURN, sensor_fusion, prev_path_size, time_inc);
              left_cost = calc_cost(car_s, car_speed*0.44704, car_lane, left_front_car, left_back_car, LTURN, NOT_TURN, sensor_fusion, prev_path_size, time_inc);
            } else if (car_lane == 0) {
              front_cost = calc_cost(car_s, car_speed*0.44704, car_lane, front_car, back_car, NOT_TURN, NOT_TURN, sensor_fusion, prev_path_size, time_inc);
              right_cost = calc_cost(car_s, car_speed*0.44704, car_lane, right_front_car, right_back_car, NOT_TURN, RTURN, sensor_fusion, prev_path_size, time_inc);
              left_cost = 99999;
            } else if (car_lane == 2) {
              front_cost = calc_cost(car_s, car_speed*0.44704, car_lane, front_car, back_car, NOT_TURN, NOT_TURN, sensor_fusion, prev_path_size, time_inc);
              right_cost = 99999;
              left_cost = calc_cost(car_s, car_speed*0.44704, car_lane, left_front_car, left_back_car, LTURN, NOT_TURN, sensor_fusion, prev_path_size, time_inc);
            }

            cout << "------------ COSTS ------------" << endl;
            cout << "this is front cost: " << front_cost << endl;
            cout << "this is left cost: " << left_cost << endl;
            cout << "This is right cost: " << right_cost << endl;
            cout << "-------------------------------" << endl;

            if (GOING_LEFT == 1) {
              if (lane == 0 && car_d >1.95 && car_d < 2.05) {
                GOING_LEFT = 0;
              } else if (lane == 1 && car_d >5.95 && car_d < 6.05) {
                GOING_LEFT = 0;
              }
            }

            if (GOING_RIGHT == 1) {
              if (lane == 2 && car_d >9.95 && car_d < 10.05) {
                GOING_RIGHT = 0;
              } else if (lane == 1 && car_d >5.95 && car_d < 6.05) {
                GOING_RIGHT = 0;
              }
            }

            if ( (front_cost < left_cost) && (front_cost < right_cost) && (GOING_LEFT == 0) && (GOING_RIGHT ==0)) {
              TURN_COOLDOWN -= 1;
              lane = car_lane;
              update_velocity(car_s, front_car, prev_path_size, time_inc);

            } else if ( (left_cost < front_cost) && (left_cost < right_cost) && (GOING_RIGHT==0) && (v_targ > 14.0)) {   // change lanes only if faster than 30 mph
              TURN_COOLDOWN = 5.0;
              lane = car_lane - 1;
              GOING_LEFT = 1;
              if (lane < 0) {
                lane = 0;
              }
              update_velocity(car_s, left_front_car, prev_path_size, time_inc);

            } else if ( (right_cost < front_cost) && (right_cost < left_cost) && (GOING_LEFT==0) && (v_targ > 14.0) ) {  // change lanes only if faster than 30 mph
              TURN_COOLDOWN = 5.0;
              lane = car_lane +1;
              GOING_RIGHT = 1;
              if (lane > 2) {
                lane = 2;
              }
              update_velocity(car_s, right_front_car, prev_path_size, time_inc);

            } else if (GOING_LEFT==0 && GOING_RIGHT==0) {
              TURN_COOLDOWN -= 1;
              lane = car_lane;
              update_velocity(car_s, front_car, prev_path_size, time_inc);
            }

            if (prev_path_size < 2) {
              // we find out what the initial car position is based on the current angle
              double initial_car_x = car_x - cos(car_yaw);
              double initial_car_y = car_y - sin(car_yaw);

              // add these points to our x and y list
              x_val.push_back(initial_car_x);
              x_val.push_back(car_x);

              y_val.push_back(initial_car_y);
              y_val.push_back(car_y);

            } else {
              start_x = previous_path_x[prev_path_size-1];
              start_y = previous_path_y[prev_path_size-1];

              double initial_start_x = previous_path_x[prev_path_size-2];
              double initial_start_y = previous_path_y[prev_path_size-2];
              start_yaw = atan2(start_y - initial_start_y, start_x - initial_start_x);

              x_val.push_back(initial_start_x);
              x_val.push_back(start_x);

              y_val.push_back(initial_start_y);
              y_val.push_back(start_y);

            }

            cout << "My current d: " << car_d << endl;
            cout << "My current lane: " << car_lane << endl;
            cout << "Target lane: " << lane << endl;
            cout << "am I going left? " << GOING_LEFT << endl;
            cout << "am i going right? " << GOING_RIGHT << endl;

            // create 3 way points 30 m away away from initial source
            // need to use lane 0 1 or 2 for left mid and right

            vector<double> waypt1 = getXY(car_s+SPACING, (2.0+4.0*(double)lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> waypt2 = getXY(car_s+SPACING*2, (2.0+4.0*(double)lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> waypt3 = getXY(car_s+SPACING*3, (2.0+4.0*(double)lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            x_val.push_back(waypt1[0]);
            x_val.push_back(waypt2[0]);
            x_val.push_back(waypt3[0]);

            y_val.push_back(waypt1[1]);
            y_val.push_back(waypt2[1]);
            y_val.push_back(waypt3[1]);

            int value_size = x_val.size();

            // convert to car coordinate system
            for(int i = 0; i<value_size; i++) {
              // find the delta between current pos(x,y) and way point(x,y)
              double dx = x_val[i] - start_x;
              double dy = y_val[i] - start_y;
              // calculate new wavepoint position relative to current point
              x_val[i] = dx*cos(start_yaw) + dy*sin(start_yaw);
              y_val[i] = dy*cos(start_yaw) - dx*sin(start_yaw);
            }

            tk::spline s;

            s.set_points(x_val, y_val);

            vector<double> next_x_vals;
            vector<double> next_y_vals;

            // initialize next x and y vals  with previous path points
            for (int i = 0; i< previous_path_x.size(); i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // break up previous path points to achieve correct velocity
            double x_targ = 30.0;
            double y_targ = s(x_targ);
            double dist_targ = sqrt((x_targ*x_targ)+(y_targ*y_targ));

            double x_inc = 0;

            // fill out the planner with new points
            for (int i = 1; i <= 50-previous_path_x.size(); i++) {

              double N = (dist_targ/(0.02*v_targ));   // speed limit is in m/s
              double x_point = x_inc + (x_targ)/N;
              double y_point = s(x_point);

              x_inc = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

              x_point = (x_ref * cos(start_yaw)-y_ref*sin(start_yaw));
              y_point = (x_ref*sin(start_yaw)+y_ref*cos(start_yaw));

              x_point += start_x;
              y_point += start_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }

            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
