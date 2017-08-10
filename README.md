# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

## Purpose
The purpose of this project is to create a planner for a simulated car to navigate high way conditions with random cars driving in random locations at random speeds.  This is implemented using C++.

## Strategy
This project was conducted in many steps, trials and errors.  It was and is a very difficult project for me to complete.  Here is the story of my discovering on how I went about doing this project.  I went about two different routes in completing this project.  Initially prior to being receiving the walkthrough video via youtube.  I spent two weeks getting my JMT version functioning.  I basically fed into this JMT the start_s, start_velocity, start_acceleration and also the end_s, end_velocity, and end_acceleration.  I also gave it a time interval T which I initially set to about 2 and tweaked up up and down.  This was repeated for also the d parameter.  

```
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
```
Then the this equation was calculated and the s and d value were calculated for each time increment of 0.02.  Then these were converted into XY coordinates and added to the list of points for the car to travel to and about.  And it was almost working but not very robust.  Also it was not able to loop the tracks correctly.  What I mean is that after about 6945, the s values restart at 0, and this for some reason was having trouble doing that.

```
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
```

Also by using the above method, somehow there would be some non smooth velocity and acceleration spikes in random parts of the track.  So, for some reason it just wasn't smooth.

Anyhow, before I even got to implementing JMT, I first downloaded and installed the simulator and got the workspace working.  Then I figured out how to first move the car down the road by simply increasing the s values while keeping the d values fixed.  I decided to start off implementing the JMT and PTG functions as shown above.

My initial approach was too complicated.  I initially tried to randomly choose different spots and calculate the cost of each trajectory and have it choose the best cost.  This resulted in having the car teleport around in random directions.  Then I figured i should just simplify this and just try to have the car move straight using the JMT function you see above.  But even this proved challenging, it wasn't moving smoothly.  

Then, after reading forums, I was able to find a way to have the car move straight in one lane using the JMT method.  Then I decided to see how to make this car change lanes.  So, i built cost functions around simply just choosing 3 different trajectories one for left, right and middle (I would later use these same cost functions when I changed the trajectory generation).

I created a speed cost, in such a way that if a car in a neighboring lane is moving faster then my car would more prefer to be in that lane.  I also created a distance cost where if there was more space in the front in a neighboring lane, then the car would prefer to be in that lane.  I also created a "bubble cost" where if there were some cars both in front and behind my car in a neighboring lane are within a certain BUBBLE_DIST of 20 meters, then the car will extremely prefer NOT to be in those lanes.  I also added a small lane change cost so that the car does not enjoy driving like a daredevil and changing lanes all the time.

These costs are calculated in the cost functions shown below:

```
// here is the lane change cost
if ((left_move ==1) || (right_move==1)) {
  total_cost += LANE_WEIGHT;
}
//cout << "total cost after lane weight: " << total_cost << endl;

// here is the speed limit costs
if (front_car_speed > SPEED_LIMIT) {
  total_cost += 0;
} else {
  total_cost += SPEED_WEIGHT*((SPEED_LIMIT-front_car_speed)/SPEED_LIMIT);
}
//cout << "total cost after speed weight: " << total_cost << endl;

// here is the distance costs
double FAR_SPACE = 120.0;
if (front_car_dist > FAR_SPACE) {
  total_cost += 0;
} else if (front_car_dist <= FAR_SPACE && front_car_dist > 0){
  total_cost += DISTANCE_WEIGHT*((FAR_SPACE-front_car_dist)/FAR_SPACE);
} else {
  total_cost += 0;
}

// Here is my collision protection costs
// this kicks in only if a left or right move is detected
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

```


Anyway, this JMT method sometimes make it all the way to the end.  But sometimes, it won't make it all the way because the path wasn't smooth.  So, now i spent many hours and days tweaking speed and time frame (parameter T for the JMT function) and play around with the cost functions to see if this would allow my car friend to make it to the goal.  And of course none of these worked very robustly.  There would also be some stretch that would cause the car to fail.  Not to mention that sometimes some random car will swurve in front of my car and brake and there's almost no way i could recover from that.  

I also found that I would hit a brick wall near the max_s = 6945.554; it will wrap back to zero.  I implemented a modulus to try to fix this:

```
/*
double wrap_distance(double s_input) {
    return s_input - max_s*floor(s_input/max_s);
} */
```

But it still didn't work.  So, then the strategy was to take a break from this project after I heard udacity was going to give out some cool pointers!  So, I took some days break until udacity released a walkthrough video which described a totally different method which I used.

So, in this method proposed by udacity.  The core idea is we want to change the global coordinates to vehicle local coordinates.  And instead of using the JMT, they use this spline tool to generate the trajectory.  So, here is the example of transforming to local vehicle coordinates:

```
// convert to car coordinate system
for(int i = 0; i<value_size; i++) {
  // find the delta between current pos(x,y) and way point(x,y)
  double dx = x_val[i] - start_x;
  double dy = y_val[i] - start_y;
  // calculate new wavepoint position relative to current point
  x_val[i] = dx*cos(start_yaw) + dy*sin(start_yaw);
  y_val[i] = dy*cos(start_yaw) - dx*sin(start_yaw);
}
```

We then create three way points about 30m apart (SPACING=30).  Note that here lane 0 is left lane, lane 1 is middle lane and lane 2 is the right lane.  And depending on the lane number the car will go to that lane in terms of d position coordinates:

```
vector<double> waypt1 = getXY(car_s+SPACING, (2.0+4.0*(double)lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
vector<double> waypt2 = getXY(car_s+SPACING*2, (2.0+4.0*(double)lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
vector<double> waypt3 = getXY(car_s+SPACING*3, (2.0+4.0*(double)lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

```

We then make use of their favorite spline tool and add points in between these way points spaced perfectly apart in a way not to violate speed limit.  In this manner, we do not have to have a speed or acceleration cost functions since the speed and acceleration is limited by how we have already spaced out these points:

```
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
```

Note that near the end shift it back to global coordinates before pushing back to the next_x_vals and next_y_vals lists.  This method provided more robust for the simulator and gave pretty smooth driving.  Therefore, I took this method and implemented my cost functions as described above to this method and the car is able to not only finish the track by also even reached up to 19.56 miles since it's running while I am writing this file.

The other neat trick was how we would slow the car down and start the car up slowly.  This was implemented in the update_velocity function:

```
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
```

Here we create this cool variable called warning where if we are too close to a vehicle then we will decrease the target velocity.  If we are too far then increase the target velocity until the speed limit was hit.  The units here is in meters per second.  Also note that by initializing the  v_targ to 0, we fix also the issue of accelerating from 0 to 47mph without any jerk.

## Conclusion

The car is able to drive around 47mph max since I made that the limit so we have some buffer zone before going over the speed limit.  Also, we should be promoting safe driving so that if there are any kids watching this they won't learn to drive fast.  :)

The car is also able to change lanes to a more empty lane.  It also stays within the lane bounds.  There are also no velocity/acceleration/jerk issues.  It is also able to go more than 4.3 miles without any incident, in fact I have seen it go up to 19.56 miles.  Therefore I consider this project a success and a very good experience.  I really learned a lot doing this project and would like to also give a special thanks to udacity for the guidance to help me finish this project.


## Future Work and possible incidents

In some random or rare occasion, possible issues would arise if a car cut my car off and braked.  Or if another car swurves into my lane while my car is there.  Or if I change lanes at the same time another car changes lanes.  I think fixing these issues 100% perfect is beyond the scope of this course.  But it should be one that should be thoroughly investigated especially if using this code in a real car.  
