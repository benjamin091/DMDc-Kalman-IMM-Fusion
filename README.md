Predictive Control of Electric Vehicles Using DMD, Kalman Filter, and IMM
Overview

Accurate prediction of power demand and vehicle speed is essential for model predictive controllers, especially in electric vehicles where recharging options are more limited compared to conventional cars.

To improve prediction quality, three different models are generated to account for the varying driving conditions, such as city driving or highway driving.

The prediction of power and speed is based on a data-driven approach using Dynamic Mode Decomposition (DMD).
Methods

Dynamic Mode Decomposition (DMD)

    DMD is a data-driven algorithm similar to matrix factorization and Principal Component Analysis (PCA).

    It processes multivariate time-series data to compute a set of dynamic modes, each with a fixed oscillation frequency and decay/growth rate.

    DMD enables interpretation of the temporal behavior of data in physically meaningful modes.

    A key feature is its ability to perform multivariate time-series prediction.

    The main idea is to find a mapping matrix AA such that X′=AXX′=AX, where XX and X′X′ represent shifted data snapshots.

    Regression steps and Singular Value Decomposition (SVD) are used to reduce noise and dimensionality.

    The resulting discrete-time DMD model approximates the system dynamics:
    x~k+1=Ax~k
    x~k+1​=Ax~k​

    The model is extended to include process and measurement noise, making it suitable for real-world noisy data.

Kalman Filter

    To improve state estimation, the DMD model is combined with a Kalman filter.

    The Kalman filter uses process and measurement noise covariances to compute an optimal gain KK.

    Discrete-time Kalman filter equations:
    x^k+1=Ax^k+K(y~k−CKalx^k)
    x^k+1​=Ax^k​+K(y~​k​−CKal​x^k​)
    y^k=CKalx^k+DKaluKal,k
    y^​k​=CKal​x^k​+DKal​uKal,k​

    Here, y~ky~​k​ are noisy measurements from DMD used as Kalman filter inputs.

Interacting Multiple Model (IMM) Approach

    Different driving conditions (city, rural, highway) require different models.

    The IMM approach manages multiple models simultaneously and switches or mixes between them to improve prediction accuracy.

    Transition and Markov matrices govern the switching behavior among models.

DMD with Control Inputs (DMDc)

    Extends DMD by incorporating control inputs via an input matrix BB.

    The system dynamics become:
    X′=AX+BY
    X′=AX+BY

    Here, YY includes inputs like position and energy change.

    The model is estimated via SVD and regression methods.

Simulation Scenarios

    Prediction quality of the IMM approach tested on transitions among three road types: "WienCity", "Landstraße", and "Autobahn".

    Parameters include transition matrix μ^μ^​, Markov matrix pp, number of regressors, prediction horizon, and noise covariances.

    Comparative studies between IMM and Kalman filter prediction accuracy under various scenarios and datasets.

Variables and Dimensions (example)
Symbol	Description	Dimensions
AA	System matrix	3j×3j3j×3j or 2j×2j2j×2j for DMDc
BB	Input matrix (DMDc only)	3j×33j×3 or 2j×22j×2
CC	Output matrix	3×3j3×3j or 2×2j2×2j
x~kx~k​	State vector after DMD estimate	3j×13j×1 or 2j×12j×1
x^kx^k​	Kalman filter state estimate	3j×13j×1 or 2j×12j×1
KK	Kalman gain	3j×33j×3 or 2j×22j×2
