# MPC Project
## The Model

I use the model that includes 6 state variables:
* x - x position of car in car coordinates (x is forward/back)
* y - y position of car in car coordinates (y is left/right)
* psi - angle of the car in car coordinates (0 is straight ahead)
* v - velocity of the car
* cte - cross-track error of the car, a.k.a. how far away from the desired center line is it.
* epsi - psi error, a.k.a. how far away is the angle from the desired angle

It also includes two actuator variables:
* delta - the steering angle to set the car to
* a - the acceleration (throttle) of the car

The update equations are as follows:

![](https://udacity-reviews-uploads.s3.amazonaws.com/_attachments/42143/1497468933/Screenshot_2017-06-14_21-35-17.png)

    x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
    y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
    psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
    v_[t+1] = v[t] + a[t] * dt
    cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
    epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt

## Timestep Length and Elapsed Duration (N & dt)

I choose to set N and dt to 5 and 0.1, respectively. This was empirically determined by manually changing up and down gradually while watching how far ahead the car was planning while seeing if computation slowed down. I chose values that cause the car to see half a second into the future, enough to understand upcoming turns at high speeds.

Higher values of dt were found to be too large to effectively steer the car, while smaller values caused oscillations.

## Polynomial Fitting and MPC Preprocessing

I preprocess speed to convert it to m/s, and waypoints by converting them into car coordinates.

I preprocess vehicle state by observing the effect of latency and advancing the vehicle state to that time in the future. This new state is used as the t0 state for the MPC.

## Model Predictive Control with Latency

    double psides0 = atan(coeffs[1]);
    double latency_dt = 1/1000.0 * ((double) config["latency_ms"]);
    double x_latency = 0 + v * cos(0) * latency_dt;
    double y_latency = 0 + v * sin(0) * latency_dt;
    double psi_latency = 0 - v / mpc.Lf * delta * latency_dt;
    double v_latency = v + a * latency_dt;
    double epsi_latency = 0 - psides0 + v * delta / mpc.Lf * latency_dt;
    double cte_latency = polyeval(coeffs, x_latency) - y_latency;