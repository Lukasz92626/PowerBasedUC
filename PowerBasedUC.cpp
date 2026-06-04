#include <ilcplex/ilocplex.h>

#include <vector>
#include <iostream>
#include <iomanip>

ILOSTLBEGIN

// Structure storing all input data for the thermal unit commitment problem
struct ThermalFleet
{
    int horizon;    // Number of time periods
    int unitCount;  // Number of generators

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

    // Minimum up/down
    std::vector<int> minOnline;             // Minimum ON time
    std::vector<int> minOffline;            // Minimum OFF time

    // Initial conditions
    std::vector<int> initialCommitment;     // Initial on/off status
    std::vector<double> initialProduction;  // Initial generation

    // System demand
    std::vector<double> load;

    // Reserve requirements
    std::vector<double> reserveUp;          // Upward reserve requirement
    std::vector<double> reserveDown;        // Downward reserve requirement
};

// Creates test data instance
ThermalFleet buildFleet()
{
    ThermalFleet f;

    f.horizon = 6;
    f.unitCount = 3;

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

    f.minOnline = { 2, 2, 1 };
    f.minOffline = { 2, 1, 1 };

    f.initialCommitment = { 0, 0, 0 };
    f.initialProduction = { 0, 0, 0 };

    f.load =
    {
        180,
        420,
        520,
        470,
        260,
        150
    };

    f.reserveUp =
    {
        40,
        50,
        60,
        50,
        30,
        20
    };

    f.reserveDown =
    {
        40,
        50,
        60,
        50,
        30,
        20
    };

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
        std::vector<std::vector<IloNumVar>> reservePositive(G, std::vector<IloNumVar>(T));
        std::vector<std::vector<IloNumVar>> reserveNegative(G, std::vector<IloNumVar>(T));

        // Variable initialization
        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                online[g][t] = IloBoolVar(env);
                startup[g][t] = IloBoolVar(env);
                shutdown[g][t] = IloBoolVar(env);

                // Power above minimum output
                aboveMinimum[g][t] = IloNumVar(
                        env,
                        0.0,
                        fleet.maximumOutput[g] - fleet.minimumOutput[g],
                        ILOFLOAT
                    );

                // Total power output
                totalOutput[g][t] = IloNumVar(
                        env,
                        0.0,
                        fleet.maximumOutput[g],
                        ILOFLOAT
                    );

                // Produced energy (integrated power)
                producedEnergy[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);

                // Reserve variables
                reservePositive[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
                reserveNegative[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            }
        }

        // -------------------------------------------------
        // (16) OBJECTIVE
        // -------------------------------------------------
        // Minimize total operating cost

        IloExpr objective(env);

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                // Variable production cost
                objective += fleet.marginalCost[g] * producedEnergy[g][t];

                // Fixed ON cost
                objective += fleet.fixedCost[g] * online[g][t];

                // Effective startup/shutdown costs (including fixed cost)
                double CSU_eff = fleet.startupCost[g] + fleet.fixedCost[g] * 1;
                double CSD_eff = fleet.shutdownCost[g] + fleet.fixedCost[g] * 1;

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
            for (int t = 1; t < T; ++t)
            {
                model.add(online[g][t] - online[g][t - 1]
                    ==
                    startup[g][t] - shutdown[g][t]
                );
            }

            // Initial condition
            model.add(online[g][0] - fleet.initialCommitment[g]
                ==
                startup[g][0] - shutdown[g][0]
            );
        }

        // -------------------------------------------------
        // (8) MINIMUM UP TIME
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = fleet.minOnline[g] - 1; t < T; ++t)
            {
                IloExpr lhs(env);

                for (int i = t - fleet.minOnline[g] + 1; i <= t; ++i)
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
            for (int t = fleet.minOffline[g] - 1; t < T; ++t)
            {
                IloExpr lhs(env);

                for (int i = t - fleet.minOffline[g] + 1; i <= t; ++i)
                {
                    lhs += shutdown[g][i];
                }

                model.add(lhs <= 1 - online[g][t]);

                lhs.end();
            }
        }

        // -------------------------------------------------
        // (6) GENERATION LIMITS
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                model.add(
                    aboveMinimum[g][t]
                    <=
                    (fleet.maximumOutput[g] - fleet.minimumOutput[g]) * online[g][t]
                );
            }
        }

        // -------------------------------------------------
        // (14) TOTAL POWER
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                model.add(totalOutput[g][t]
                    ==
                    fleet.minimumOutput[g] * online[g][t] + aboveMinimum[g][t]
                );
            }
        }

        // -------------------------------------------------
        // (15) ENERGY EQUATIONS
        // -------------------------------------------------
        // Trapezoidal approximation

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                IloExpr avg(env);

                avg += totalOutput[g][t];

                if (t > 0)
                {
                    avg += totalOutput[g][t - 1];
                }
                else
                {
                    avg += fleet.initialProduction[g];
                }

                model.add(producedEnergy[g][t] == 0.5 * avg);

                avg.end();
            }
        }

        // -------------------------------------------------
        // (27) RAMP-UP
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 1; t < T; ++t)
            {
                model.add(totalOutput[g][t] - totalOutput[g][t - 1]
                    <=
                    fleet.rampUpLimit[g] * online[g][t - 1] + fleet.startupCap[g] * startup[g][t]
                );
            }
        }

        for (int g = 0; g < G; ++g)
        {
            model.add(
                totalOutput[g][0] - fleet.initialProduction[g]
                <=
                fleet.rampUpLimit[g] * fleet.initialCommitment[g] + fleet.startupCap[g] * startup[g][0]
            );
        }

        // -------------------------------------------------
        // (28) RAMP-DOWN
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 1; t < T; ++t)
            {
                model.add(totalOutput[g][t - 1]  - totalOutput[g][t]
                    <=
                    fleet.rampDownLimit[g] * online[g][t] + fleet.shutdownCap[g] * shutdown[g][t]
                );
            }
        }

        for (int g = 0; g < G; ++g)
        {
            model.add(
                fleet.initialProduction[g] - totalOutput[g][0]
                <=
                fleet.rampDownLimit[g] * online[g][0] + fleet.shutdownCap[g] * shutdown[g][0]
            );
        }

        // -------------------------------------------------
        // (AUX) RESERVE RAMP LIMIT (UP)
        // -------------------------------------------------
        // Reserve cannot exceed ramping capability

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                model.add(
                    reservePositive[g][t] <= fleet.rampUpLimit[g] * online[g][t]
                );
            }
        }

        // -------------------------------------------------
        // (26) UP RESERVE CAPABILITY
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T - 1; ++t)
            {
                model.add(
                    aboveMinimum[g][t] + reservePositive[g][t]
                    <=
                    (fleet.maximumOutput[g] - fleet.minimumOutput[g]) * online[g][t]
                    - (fleet.maximumOutput[g] - fleet.shutdownCap[g]) * shutdown[g][t + 1]
                    + (fleet.startupCap[g] - fleet.minimumOutput[g]) * startup[g][t + 1]
                );
            }
        }

        int t_last = T - 1;

        for (int g = 0; g < G; ++g)
        {
            model.add(
                aboveMinimum[g][t_last] + reservePositive[g][t_last]
                <=
                (fleet.maximumOutput[g] - fleet.minimumOutput[g]) * online[g][t_last]
            );
        }

        // -------------------------------------------------
        // (27) (28) DOWN RESERVE CAPABILITY
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                // (28)
                model.add(
                    reserveNegative[g][t]
                    <=
                    totalOutput[g][t] - fleet.minimumOutput[g] * online[g][t] -
                    ( fleet.shutdownCap[g] - fleet.minimumOutput[g]) * shutdown[g][t]
                );

                // (27)
                model.add(
                    reserveNegative[g][t]
                    <=
                    totalOutput[g][t] - fleet.minimumOutput[g] * online[g][t]
                );
            }
        }

        // -------------------------------------------------
        // (19) SYSTEM DEMAND BALANCE
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
            IloExpr reserve(env);

            for (int g = 0; g < G; ++g)
            {
                reserve += reservePositive[g][t];
            }

            model.add(reserve >= fleet.reserveUp[t]);

            reserve.end();
        }

        // -------------------------------------------------
        // (21) DOWN RESERVE REQUIREMENT
        // -------------------------------------------------

        for (int t = 0; t < T; ++t)
        {
            IloExpr reserve(env);

            for (int g = 0; g < G; ++g)
            {
                reserve += reserveNegative[g][t];
            }

            model.add(reserve >= fleet.reserveDown[t]);

            reserve.end();
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

            std::cout
                << "t"
                << "\tON"
                << "\tSU"
                << "\tSD"
                << "\tP"
                << "\tRU"
                << "\tRD"
                << "\tE\n";

            for (int t = 0; t < T; ++t)
            {
                std::cout
                    << t
                    << "\t"
                    << solver.getValue(online[g][t])
                    << "\t"
                    << solver.getValue(startup[g][t])
                    << "\t"
                    << solver.getValue(shutdown[g][t])
                    << "\t"
                    << solver.getValue(totalOutput[g][t])
                    << "\t"
                    << solver.getValue(reservePositive[g][t])
                    << "\t"
                    << solver.getValue(reserveNegative[g][t])
                    << "\t"
                    << solver.getValue(producedEnergy[g][t])
                    << "\n";
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