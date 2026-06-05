#include <ilcplex/ilocplex.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

ILOSTLBEGIN

// Structure storing all input data for the thermal unit commitment problem
struct ThermalFleet
{
    int horizon;    // Number of time periods
    int unitCount;  // Number of generators

    std::vector<bool> isQuickStart;

    // Costs parameters
    std::vector<double> marginalCost;       // Variable production cost [$/MWh]
    std::vector<double> fixedCost;          // Fixed on/off cost [$/h]
    std::vector<double> startupCost;        // Startup cost [$]
    std::vector<double> shutdownCost;       // Shutdown cost [$]

    // Power limits
    std::vector<double> minimumOutput;      // Minimum generation level
    std::vector<double> maximumOutput;      // Maximum generation level

    // Ramping limits
    std::vector<double> rampUpLimit;        // Maximum ramp-up per period
    std::vector<double> rampDownLimit;      // Maximum ramp-down per period

    // Startup / shutdown capability
    std::vector<double> startupCap;         // Max output during startup
    std::vector<double> shutdownCap;        // Max output during shutdown

    // Slow-start: duration of startup/shutdown process
    std::vector<int> startupDuration;       // SUD_g [h]
    std::vector<int> shutdownDuration;      // SDD_g [h]

    // Slow-start: power trajectory during startup PSU_i
    // Outer vector: generator g; inner vector: SUD_g values
    std::vector<std::vector<double>> startupTraj;   // PSU_{g,i}

    // Slow-start: power trajectory during shutdown PSD_i
    std::vector<std::vector<double>> shutdownTraj; // PSD_{g,i}

    // Minimum up/down
    std::vector<int> minOnline;             // Minimum ON time
    std::vector<int> minOffline;            // Minimum OFF time

    // Initial conditions
    std::vector<int> initialStatus;         // Initial on/off status
    std::vector<double> initialProduction;  // Initial generation

    // System demand
    std::vector<double> load;               // D_t [MW]

    // Reserve requirements
    std::vector<double> reserveUp;          // Upward reserve requirement
    std::vector<double> reserveDown;        // Downward reserve requirement
};

// Build a linear startup power trajectory for slow-start units.
static std::vector<double> makeStartupTraj(double Pmin, int SUD)
{
    std::vector<double> traj(SUD);
    for (int i = 0; i < SUD; ++i)
        traj[i] = Pmin * (i + 1) / static_cast<double>(SUD + 1);
    return traj;
}

// Build a linear shutdown power trajectory for slow-start units.
static std::vector<double> makeShutdownTraj(double Pmin, int SDD)
{
    // PSD_i  for i = 1 … SDD+1  where PSD_1 = P̲, PSD_{SDD+1} = 0
    // We store i = 1 … SDD+1 (size SDD+1) so indexing matches eq.(12).
    std::vector<double> traj(SDD + 1);
    for (int i = 0; i <= SDD; ++i)
        traj[i] = Pmin * (SDD - i) / static_cast<double>(SDD);
    return traj;
}

// Creates test data instance
ThermalFleet buildFleet()
{
    ThermalFleet f;

    f.horizon = 6;
    f.unitCount = 3;

    f.isQuickStart = { true, true, false };

    f.marginalCost = { 18, 24, 40 };
    f.fixedCost = { 120, 80, 30 };
    f.startupCost = { 600, 400, 100 };
    f.shutdownCost = { 150, 120, 30 };

    f.minimumOutput = { 120, 80, 40 };
    f.maximumOutput = { 350, 200, 120 };

    f.rampUpLimit = { 120, 90, 80 };
    f.rampDownLimit = { 120, 90, 80 };

    f.startupCap = { 220, 150, 100 };
    f.shutdownCap = { 220, 150, 100 };

    f.startupDuration = { 1, 1, 2 };
    f.shutdownDuration = { 1, 1, 2 };

    f.startupTraj = { {}, {}, makeStartupTraj(40.0, 2) };
    f.shutdownTraj = { {}, {}, makeShutdownTraj(40.0, 2) };

    f.minOnline = { 2, 2, 1 };
    f.minOffline = { 2, 1, 1 };

    f.initialStatus = { 0, 0, 0 };
    f.initialProduction = { 0, 0, 0 };

    f.load = { 180, 420, 520, 470, 260, 150 };

    f.reserveUp = { 40, 50, 60, 50, 30, 20 };
    f.reserveDown = { 40, 50, 60, 50, 30, 20 };

    return f;
}

int main()
{
    IloEnv env;

    try
    {
        ThermalFleet fleet = buildFleet();

        IloModel model(env);

        const int G = fleet.unitCount;  // number of generators
        const int T = fleet.horizon;    // number of time periods

        // -------------------------------------------------
        // VARIABLES
        // -------------------------------------------------

        // Binary commitment variables
        std::vector<std::vector<IloBoolVar>> online(G, std::vector<IloBoolVar>(T));
        std::vector<std::vector<IloBoolVar>> startup(G, std::vector<IloBoolVar>(T));
        std::vector<std::vector<IloBoolVar>> shutdown(G, std::vector<IloBoolVar>(T));

        // Continuous variables
        std::vector<std::vector<IloNumVar>> aboveMinimum(G, std::vector<IloNumVar>(T));
        std::vector<std::vector<IloNumVar>> totalOutput(G, std::vector<IloNumVar>(T));
        std::vector<std::vector<IloNumVar>> producedEnergy(G, std::vector<IloNumVar>(T));
        std::vector<std::vector<IloNumVar>> reservePos(G, std::vector<IloNumVar>(T));
        std::vector<std::vector<IloNumVar>> reserveNeg(G, std::vector<IloNumVar>(T));

        // Variable initialization
        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                online[g][t] = IloBoolVar(env);
                startup[g][t] = IloBoolVar(env);
                shutdown[g][t] = IloBoolVar(env);

                // Power above minimum output
                aboveMinimum[g][t] = IloNumVar(env, 0.0, fleet.maximumOutput[g] - fleet.minimumOutput[g], ILOFLOAT);
                
                // Total power output
                totalOutput[g][t] = IloNumVar(env, 0.0, fleet.maximumOutput[g], ILOFLOAT);

                // Produced energy (integrated power)
                producedEnergy[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);

                // Reserve variables
                reservePos[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
                reserveNeg[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            }
        }

        // -------------------------------------------------
        // (16) OBJECTIVE
        // -------------------------------------------------
        // Minimize total operating cost

        IloExpr objective(env);

        for (int g = 0; g < G; ++g)
        {
            int SUD = fleet.startupDuration[g];
            int SDD = fleet.shutdownDuration[g];

            // Effective startup/shutdown costs (including fixed cost)
            double CSU_eff = fleet.startupCost[g] + fleet.fixedCost[g] * SUD;
            double CSD_eff = fleet.shutdownCost[g] + fleet.fixedCost[g] * SDD;

            for (int t = 0; t < T; ++t)
            {
                // Variable production cost
                objective += fleet.marginalCost[g] * producedEnergy[g][t];

                // Fixed ON cost
                objective += fleet.fixedCost[g] * online[g][t];

                objective += CSU_eff * startup[g][t];
                objective += CSD_eff * shutdown[g][t];
            }
        }

        model.add(IloMinimize(env, objective));

        // -------------------------------------------------
        // (7) LOGIC EQUATIONS
        // -------------------------------------------------
        // Commitment state transitions

        for (int g = 0; g < G; ++g)
        {
            // t = 0: use initialStatus to determine u_{-1}
            int u0_prev = (fleet.initialStatus[g] > 0) ? 1 : 0;
            model.add(online[g][0] - u0_prev == startup[g][0] - shutdown[g][0]);

            for (int t = 1; t < T; ++t) {
                model.add(online[g][t] - online[g][t - 1] == startup[g][t] - shutdown[g][t]);
            }
        }

        // -------------------------------------------------
        // (8) MINIMUM UP TIME
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            int TU = fleet.minOnline[g];

            // Initial commitment constraint for units already online
            if (fleet.initialStatus[g] > 0)
            {
                int hoursOn = fleet.initialStatus[g];
                int mustRemainOn = std::max(0, TU - hoursOn);

                for (int t = 0; t < std::min(mustRemainOn, T); ++t) {
                    model.add(online[g][t] == 1);
                }
            }

            // Rolling minimum-uptime constraint
            for (int t = TU - 1; t < T; ++t)
            {
                IloExpr lhs(env);

                for (int i = t - TU + 1; i <= t; ++i)
                {
                    lhs += startup[g][i];
                }

                model.add(lhs <= online[g][t]);
                lhs.end();
            }
        }

        // -------------------------------------------------
        // (9) MINIMUM DOWN TIME
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            int TD = fleet.minOffline[g];

            // Initial commitment constraint for units already offline
            if (fleet.initialStatus[g] < 0)
            {
                int hoursOff = -fleet.initialStatus[g];
                int mustRemainOff = std::max(0, TD - hoursOff);

                for (int t = 0; t < std::min(mustRemainOff, T); ++t)
                    model.add(online[g][t] == 0);
            }

            // Rolling minimum-downtime constraint
            for (int t = TD - 1; t < T; ++t)
            {
                IloExpr lhs(env);

                for (int i = t - TD + 1; i <= t; ++i)
                {
                    lhs += shutdown[g][i];
                }

                model.add(lhs <= 1 - online[g][t]);
                lhs.end();
            }
        }

        // -------------------------------------------------
        // (4) (5) (6) GENERATION LIMITS WITH SU/SD CAPABILITIES
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            double Pmax = fleet.maximumOutput[g];
            double Pmin = fleet.minimumOutput[g];
            double SU = fleet.startupCap[g];
            double SD = fleet.shutdownCap[g];

            for (int t = 0; t < T - 1; ++t)
            {
                // (4)
                model.add(
                    aboveMinimum[g][t] + reservePos[g][t]
                    <=
                    (Pmax - Pmin) * online[g][t]
                    - (Pmax - SD) * shutdown[g][t + 1]
                    + (SU - Pmin) * startup[g][t + 1]
                );
            }

            // (5)
            model.add(
                aboveMinimum[g][T - 1] + reservePos[g][T - 1]
                <=
                (Pmax - Pmin) * online[g][T - 1]
            );
        }

        // -------------------------------------------------------------------
        // DOWN-RESERVE CAPABILITY (eq. 26 lower part + shutdown correction)
        // -------------------------------------------------------------------
        for (int g = 0; g < G; ++g)
        {
            double Pmin = fleet.minimumOutput[g];
            double SD = fleet.shutdownCap[g];

            for (int t = 0; t < T; ++t)
            {
                // With shutdown correction
                model.add(
                    reserveNeg[g][t]
                    <=
                    aboveMinimum[g][t] - (SD - Pmin) * shutdown[g][t]
                );

                // General: reserve cannot exceed power above minimum
                model.add(reserveNeg[g][t] <= aboveMinimum[g][t]);
            }
        }

        // -------------------------------------------------------------------
        // (14) (29) TOTAL POWER for QUICK-START units
        // -------------------------------------------------------------------
        for (int g = 0; g < G; ++g)
        {
            double Pmin = fleet.minimumOutput[g];

            for (int t = 0; t < T; ++t)
            {
                IloExpr rhs(env);

                // Common term
                rhs += Pmin * online[g][t];
                if (t + 1 < T) {
                    rhs += Pmin * startup[g][t + 1];
                }

                // Controllable part above minimum
                rhs += aboveMinimum[g][t];

                if (!fleet.isQuickStart[g])
                {
                    int SUD = fleet.startupDuration[g];
                    int SDD = fleet.shutdownDuration[g];

                    // Startup trajectory term (iii) in eq. (12):
                    for (int i = 1; i <= SUD; ++i)
                    {
                        int vIdx = t - i + SUD + 2 - 1; // convert to 0-based
                        double PSU_i = fleet.startupTraj[g][i - 1];
                        if (vIdx >= 0 && vIdx < T)
                            rhs += PSU_i * startup[g][vIdx];
                    }

                    // Shutdown trajectory term (ii) in eq. (12):
                    for (int i = 2; i <= SDD + 1; ++i)
                    {
                        int wIdx = t - i + 2 - 1; // convert to 0-based
                        double PSD_i = fleet.shutdownTraj[g][i - 1]; // PSD_i
                        if (wIdx >= 0 && wIdx < T)
                            rhs += PSD_i * shutdown[g][wIdx];
                    }
                }

                model.add(totalOutput[g][t] == rhs);
                rhs.end();
            }
        }

        // -------------------------------------------------------------------
        // (15) (31) ENERGY — trapezoidal rule
        // -------------------------------------------------------------------
        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                IloExpr avg(env);
                avg += totalOutput[g][t];

                if (t > 0)
                    avg += totalOutput[g][t - 1];
                else
                    avg += fleet.initialProduction[g];

                model.add(producedEnergy[g][t] == 0.5 * avg);
                avg.end();
            }
        }

        // -------------------------------------------------------------------
        // (27) RAMP-UP WITH RESERVES
        // -------------------------------------------------------------------
        for (int g = 0; g < G; ++g)
        {
            double Pmin = fleet.minimumOutput[g];
            double RU = fleet.rampUpLimit[g];
            double SU = fleet.startupCap[g];

            // previous aboveMinimum derived from initial production
            int u0_prev = (fleet.initialStatus[g] > 0) ? 1 : 0;
            double p0_prev = std::max(0.0, fleet.initialProduction[g] - Pmin * u0_prev);

            model.add(
                aboveMinimum[g][0] + reservePos[g][0] - p0_prev
                <=
                RU * u0_prev + (SU - Pmin) * startup[g][0]
            );

            for (int t = 1; t < T; ++t)
            {
                model.add(
                    aboveMinimum[g][t] + reservePos[g][t] - aboveMinimum[g][t - 1]
                    <=
                    RU * online[g][t - 1] + (SU - Pmin) * startup[g][t]
                );
            }
        }

        // -------------------------------------------------------------------
        // (28) RAMP-DOWN WITH RESERVES
        // -------------------------------------------------------------------
        for (int g = 0; g < G; ++g)
        {
            double Pmin = fleet.minimumOutput[g];
            double RD = fleet.rampDownLimit[g];
            double SD = fleet.shutdownCap[g];

            int u0_prev = (fleet.initialStatus[g] > 0) ? 1 : 0;
            double p0_prev = std::max(0.0, fleet.initialProduction[g] - Pmin * u0_prev);

            model.add(
                p0_prev - (aboveMinimum[g][0] - reserveNeg[g][0])
                <=
                RD * online[g][0] + (SD - Pmin) * shutdown[g][0]
            );

            for (int t = 1; t < T; ++t)
            {
                model.add(
                    aboveMinimum[g][t - 1] - (aboveMinimum[g][t] - reserveNeg[g][t])
                    <=
                    RD * online[g][t] + (SD - Pmin) * shutdown[g][t]
                );
            }
        }

        // -------------------------------------------------
        // (19) SYSTEM POWER DEMAND BALANCE
        // -------------------------------------------------

        for (int t = 0; t < T; ++t)
        {
            IloExpr generation(env);

            for (int g = 0; g < G; ++g)
            {
                generation += totalOutput[g][t];
            }

            model.add(generation >= fleet.load[t]);
            generation.end();
        }

        // -------------------------------------------------
        // (20) UP RESERVE REQUIREMENT
        // -------------------------------------------------

        for (int t = 0; t < T; ++t)
        {
            IloExpr res(env);

            for (int g = 0; g < G; ++g)
            {
                res += reservePos[g][t];
            }

            model.add(res >= fleet.reserveUp[t]);
            res.end();
        }

        // -------------------------------------------------
        // (21) DOWN RESERVE REQUIREMENT
        // -------------------------------------------------

        for (int t = 0; t < T; ++t)
        {
            IloExpr res(env);

            for (int g = 0; g < G; ++g)
            {
                res += reserveNeg[g][t];
            }

            model.add(res >= fleet.reserveDown[t]);
            res.end();
        }

        // -------------------------------------------------
        // SOLVER
        // -------------------------------------------------

        IloCplex solver(model);

        solver.setParam(IloCplex::Param::Threads, 1);

        if (!solver.solve())
        {
            std::cout << "No feasible solution\n";

            env.end();

            return 1;
        }

        // -------------------------------------------------
        // RESULTS
        // -------------------------------------------------

        std::cout << std::fixed << std::setprecision(2);

        std::cout << "\n";
        std::cout << "Objective = " << solver.getObjValue() << "\n\n";

        for (int g = 0; g < G; ++g)
        {
            std::cout << "============================\n";
            std::cout << "Generator " << g << "\n";
            std::cout << "============================\n";

            std::cout << "t" << "\tON" << "\tSU" << "\tSD" << "\tp_abs" << "\tp_tot" << "\tr+" << "\tr-" << "\tE\n";

            for (int t = 0; t < T; ++t)
            {
                std::cout
                    << t << "\t"
                    << solver.getValue(online[g][t]) << "\t"
                    << solver.getValue(startup[g][t]) << "\t"
                    << solver.getValue(shutdown[g][t]) << "\t"
                    << solver.getValue(aboveMinimum[g][t]) << "\t"
                    << solver.getValue(totalOutput[g][t]) << "\t"
                    << solver.getValue(reservePos[g][t]) << "\t"
                    << solver.getValue(reserveNeg[g][t]) << "\t"
                    << solver.getValue(producedEnergy[g][t]) << "\n";
            }

            std::cout << std::endl;
        }

        objective.end();
    }
    catch (IloException& ex)
    {
        std::cerr << "CPLEX exception: " << ex << std::endl;
    }
    catch (...)
    {
        std::cerr << "Unknown exception" << std::endl;
    }

    env.end();

    return 0;
}