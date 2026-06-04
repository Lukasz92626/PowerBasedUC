#include <ilcplex/ilocplex.h>

#include <vector>
#include <iostream>
#include <iomanip>

ILOSTLBEGIN

struct ThermalFleet
{
    int horizon;
    int unitCount;

    // Costs
    std::vector<double> marginalCost;
    std::vector<double> fixedCost;
    std::vector<double> startupCost;
    std::vector<double> shutdownCost;

    // Power limits
    std::vector<double> minimumOutput;
    std::vector<double> maximumOutput;

    // Ramping
    std::vector<double> rampUpLimit;
    std::vector<double> rampDownLimit;

    // Startup / shutdown capability
    std::vector<double> startupCap;
    std::vector<double> shutdownCap;

    //std::vector<double> startupTrajectory;
    //std::vector<double> shutdownTrajectory;

    // Minimum up/down
    std::vector<int> minOnline;
    std::vector<int> minOffline;

    // Initial conditions
    std::vector<int> initialCommitment;
    std::vector<double> initialProduction;

    // Demand
    std::vector<double> load;

    // Reserve requirements
    std::vector<double> reserveUp;
    std::vector<double> reserveDown;
};

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

    //f.startupTrajectory = { 100, 70, 40 };
    //f.shutdownTrajectory = { 100, 60, 30 };

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

        const int G = fleet.unitCount;
        const int T = fleet.horizon;

        // -------------------------------------------------
        // VARIABLES
        // -------------------------------------------------

        std::vector<std::vector<IloBoolVar>> online(
            G,
            std::vector<IloBoolVar>(T)
        );

        std::vector<std::vector<IloBoolVar>> startup(
            G,
            std::vector<IloBoolVar>(T)
        );

        std::vector<std::vector<IloBoolVar>> shutdown(
            G,
            std::vector<IloBoolVar>(T)
        );

        std::vector<std::vector<IloNumVar>> aboveMinimum(
            G,
            std::vector<IloNumVar>(T)
        );

        std::vector<std::vector<IloNumVar>> totalOutput(
            G,
            std::vector<IloNumVar>(T)
        );

        std::vector<std::vector<IloNumVar>> producedEnergy(
            G,
            std::vector<IloNumVar>(T)
        );

        std::vector<std::vector<IloNumVar>> reservePositive(
            G,
            std::vector<IloNumVar>(T)
        );

        std::vector<std::vector<IloNumVar>> reserveNegative(
            G,
            std::vector<IloNumVar>(T)
        );



        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                online[g][t] = IloBoolVar(env);
                startup[g][t] = IloBoolVar(env);
                shutdown[g][t] = IloBoolVar(env);

                aboveMinimum[g][t] = IloNumVar(
                        env,
                        0.0,
                        fleet.maximumOutput[g] - fleet.minimumOutput[g],
                        ILOFLOAT
                    );

                totalOutput[g][t] = IloNumVar(
                        env,
                        0.0,
                        fleet.maximumOutput[g],
                        ILOFLOAT
                    );

                producedEnergy[g][t] =
                    IloNumVar(
                        env,
                        0.0,
                        IloInfinity,
                        ILOFLOAT
                    );

                reservePositive[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);

                reserveNegative[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            }
        }

        // -------------------------------------------------
        // OBJECTIVE
        // -------------------------------------------------

        IloExpr objective(env);

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                objective += fleet.marginalCost[g] * producedEnergy[g][t];

                objective += fleet.fixedCost[g] * online[g][t];

                //objective += fleet.startupCost[g] * startup[g][t];

                //objective += fleet.shutdownCost[g] * shutdown[g][t];

                double CSU_eff = fleet.startupCost[g] + fleet.fixedCost[g] * 1;
                double CSD_eff = fleet.shutdownCost[g] + fleet.fixedCost[g] * 1;

                objective += CSU_eff * startup[g][t];
                objective += CSD_eff * shutdown[g][t];
            }
        }

        model.add(IloMinimize(env, objective));

        // -------------------------------------------------
        // LOGIC EQUATIONS
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 1; t < T; ++t)
            {
                model.add(online[g][t] - online[g][t - 1]
                    ==
                    startup[g][t] - shutdown[g][t]
                );
            }

            model.add(online[g][0] - fleet.initialCommitment[g]
                ==
                startup[g][0] - shutdown[g][0]
            );
        }

        // -------------------------------------------------
        // MINIMUM UP TIME
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
        // MINIMUM DOWN TIME
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
        // GENERATION LIMITS
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
        // TOTAL POWER
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
        // ENERGY EQUATIONS
        // -------------------------------------------------

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
        // RAMP-UP
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
        // RAMP-DOWN
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
        // UP RESERVE CAPABILITY
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                model.add(
                    totalOutput[g][t] + reservePositive[g][t]
                    <=
                    fleet.maximumOutput[g] * online[g][t] -
                    (fleet.maximumOutput[g] - fleet.startupCap[g]) * startup[g][t]
                );

                model.add(
                    reservePositive[g][t]
                    <=
                    fleet.maximumOutput[g] - totalOutput[g][t]
                );
            }
        }

        // -------------------------------------------------
        // DOWN RESERVE CAPABILITY
        // -------------------------------------------------

        for (int g = 0; g < G; ++g)
        {
            for (int t = 0; t < T; ++t)
            {
                model.add(
                    reserveNegative[g][t]
                    <=
                    totalOutput[g][t] - fleet.minimumOutput[g] * online[g][t] -
                    ( fleet.shutdownCap[g] - fleet.minimumOutput[g]) * shutdown[g][t]
                );

                model.add(
                    reserveNegative[g][t]
                    <=
                    totalOutput[g][t] - fleet.minimumOutput[g] * online[g][t]
                );
            }
        }

        // -------------------------------------------------
        // SYSTEM DEMAND BALANCE
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
        // UP RESERVE REQUIREMENT
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
        // DOWN RESERVE REQUIREMENT
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
