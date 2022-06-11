// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Unscented Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.    

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;
// -------------------------------------------------- //

// ----------------------------------------------------------------------- //
// USEFUL HELPER FUNCTIONS
VectorXd normaliseState(VectorXd state)
{
    state(2) = wrapAngle(state(2));
    return state;
}
VectorXd normaliseLidarMeasurement(VectorXd meas)
{
    meas(1) = wrapAngle(meas(1));
    return meas;
}
std::vector<VectorXd> generateSigmaPoints(VectorXd state, MatrixXd cov)
{
    std::vector<VectorXd> sigmaPoints;

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE
    int n = state.size();
    int k = 3 - n;
    MatrixXd L = cov.llt().matrixL();

    sigmaPoints.emplace_back(state);
    for (int i = 0; i < n; i++) {
        sigmaPoints.emplace_back(state + sqrt(n + k) * L.col(i));
        sigmaPoints.emplace_back(state - sqrt(n + k) * L.col(i));
    }

    // ----------------------------------------------------------------------- //

    return sigmaPoints;
}

std::vector<double> generateSigmaWeights(unsigned int numStates)
{
    std::vector<double> weights;

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE

    double k = 3.0 - numStates;
    double wi = 0.5 / (numStates + k);
    weights.push_back(k / (numStates + k));
    for (int i = 0; i < 2 * numStates; i++) {
        weights.push_back(wi);
    }
    
    // ----------------------------------------------------------------------- //

    return weights;
}

VectorXd lidarMeasurementModel(VectorXd aug_state, double beaconX, double beaconY)
{
    VectorXd z_hat = VectorXd::Zero(2);

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE

    double x = aug_state(0);
    double y = aug_state(1);
    double psi = aug_state(2);
    double range_noise = aug_state(4);
    double bearing_noise = aug_state(5);

    double diff_x = beaconX - x;
    double diff_y = beaconY - y;
    double r = sqrt(diff_x*diff_x + diff_y*diff_y) + range_noise;
    double theta = atan2(diff_y, diff_x) - psi + bearing_noise;
    z_hat << r, theta;

    // ----------------------------------------------------------------------- //

    return z_hat;
}

VectorXd vehicleProcessModel(VectorXd aug_state, double psi_dot, double dt)
{
    VectorXd new_state = VectorXd::Zero(4);

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE

    double x = aug_state(0);
    double y = aug_state(1);
    double psi = aug_state(2);
    double v = aug_state(3);
    double heading_noise = aug_state(4);
    double acc_noise = aug_state(5);

    new_state(0) = x + dt * v * cos(psi);
    new_state(1) = y + dt * v * sin(psi);
    new_state(2) = psi + dt * (psi_dot + heading_noise);
    new_state(3) = v + dt * acc_noise;

    // ----------------------------------------------------------------------- //

    return new_state;
}
// ----------------------------------------------------------------------- //

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Update Step for the Lidar Measurements in the 
        // section below.
        // HINT: Use the normaliseState() and normaliseLidarMeasurement() functions
        // to always keep angle values within correct range.
        // HINT: Do not normalise during sigma point calculation!
        // HINT: You can use the constants: LIDAR_RANGE_STD, LIDAR_THETA_STD
        // HINT: The mapped-matched beacon position can be accessed by the variables
        // map_beacon.x and map_beacon.y
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        BeaconData map_beacon = map.getBeaconWithId(meas.id); // Match Beacon with built in Data Association Id
        if (meas.id != -1 && map_beacon.id != -1) // Check that we have a valid beacon match
        {

            MatrixXd R = Matrix2d::Zero();
            R(0, 0) = LIDAR_RANGE_STD * LIDAR_RANGE_STD;
            R(1, 1) = LIDAR_THETA_STD * LIDAR_THETA_STD;

            int n_v = 2;
            int n_z = 2;
            int n_x = state.size();
            int n_aug = n_x + n_v;
            VectorXd aug_state = VectorXd::Zero(n_aug);
            aug_state.head(n_x) = state;
            MatrixXd aug_cov = MatrixXd::Zero(n_aug, n_aug);
            aug_cov.topLeftCorner(n_x, n_x) = cov;
            aug_cov.bottomRightCorner(n_v, n_v) = R;

            std::vector<double> weights = generateSigmaWeights(n_aug);
            std::vector<VectorXd> sigma_points = generateSigmaPoints(aug_state, aug_cov);

            std::vector<VectorXd> z_sig;
            for (const auto& p : sigma_points) {
                z_sig.emplace_back(lidarMeasurementModel(p, map_beacon.x, map_beacon.y));
            }

            VectorXd z_mean = Vector2d::Zero(n_z);
            for (int i = 0; i < z_sig.size(); i++) {
                z_mean += weights[i] * z_sig[i];
            }

            MatrixXd S = MatrixXd::Zero(n_z, n_z);
            for (int i = 0; i < z_sig.size(); i++) {
                VectorXd diff = normaliseLidarMeasurement(z_sig[i] - z_mean);
                S += weights[i] * diff * diff.transpose();
            }

            Vector2d z{meas.range, meas.theta};
            VectorXd mu = normaliseLidarMeasurement(z - z_mean);

            MatrixXd cov_xz = MatrixXd::Zero(n_x, n_z);
            for (int i = 0; i < z_sig.size(); i++) {
                VectorXd diff_x = normaliseState(sigma_points[i].head(n_x) - state);
                VectorXd diff_z = normaliseLidarMeasurement(z_sig[i] - z_mean);
                cov_xz += weights[i] * diff_x * diff_z.transpose();
            }

            MatrixXd K = cov_xz * S.inverse();
            state = state + K * mu;
            cov = cov - K * S * K.transpose();

        }
        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Prediction Step for the system in the  
        // section below.
        // HINT: Assume the state vector has the form [PX, PY, PSI, V].
        // HINT: Use the Gyroscope measurement as an input into the prediction step.
        // HINT: You can use the constants: ACCEL_STD, GYRO_STD
        // HINT: Use the normaliseState() function to always keep angle values within correct range.
        // HINT: Do NOT normalise during sigma point calculation!
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        MatrixXd Q = Matrix2d::Zero();
        Q(0, 0) = GYRO_STD * GYRO_STD;
        Q(1, 1) = ACCEL_STD * ACCEL_STD;

        int n_w = 2;
        int n_x = state.size();
        int n_aug = n_x + n_w;
        VectorXd aug_state = VectorXd::Zero(n_aug);
        aug_state.head(n_x) = state;
        MatrixXd aug_cov = MatrixXd::Zero(n_aug, n_aug);
        aug_cov.topLeftCorner(n_x, n_x) = cov;
        aug_cov.bottomRightCorner(n_w, n_w) = Q;

        std::vector<double> weights = generateSigmaWeights(n_aug);
        std::vector<VectorXd> sigma_points = generateSigmaPoints(aug_state, aug_cov);

        std::vector<VectorXd> transformed;
        for (const auto& p : sigma_points) {
            transformed.emplace_back(vehicleProcessModel(p, gyro.psi_dot, dt));
        }

        state = Vector4d::Zero();
        for (int i = 0; i < transformed.size(); i++) {
            state += weights[i] * transformed[i];
        }
        state = normaliseState(state);

        cov = Matrix4d::Zero();
        for (int i = 0; i < transformed.size(); i++) {
            VectorXd diff = normaliseState(transformed[i] - state);
            cov += weights[i] * diff * diff.transpose();
        }
        
        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the UKF update state would just produce the same result.
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        VectorXd z = Vector2d::Zero();
        MatrixXd H = MatrixXd(2,4);
        MatrixXd R = Matrix2d::Zero();

        z << meas.x,meas.y;
        H << 1,0,0,0,0,1,0,0;
        R(0,0) = GPS_POS_STD*GPS_POS_STD;
        R(1,1) = GPS_POS_STD*GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        MatrixXd K = cov*H.transpose()*S.inverse();

        state = state + K*y;
        cov = (MatrixXd::Identity(4,4) - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        // You may modify this initialisation routine if you can think of a more
        // robust and accuracy way of initialising the filter.
        // ----------------------------------------------------------------------- //
        // YOU ARE FREE TO MODIFY THE FOLLOWING CODE HERE

        VectorXd state = Vector4d::Zero();
        MatrixXd cov = Matrix4d::Zero();

        state(0) = meas.x;
        state(1) = meas.y;
        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = GPS_POS_STD*GPS_POS_STD;
        cov(2,2) = INIT_PSI_STD*INIT_PSI_STD;
        cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;

        setState(state);
        setCovariance(cov);

        // ----------------------------------------------------------------------- //
    }             
}

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance()
{
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0){pos_cov << cov(0,0), cov(0,1), cov(1,0), cov(1,1);}
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState()
{
    if (isInitialised())
    {
        VectorXd state = getState(); // STATE VECTOR [X,Y,PSI,V,...]
        return VehicleState(state[0],state[1],state[2],state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt){}
