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

    bool selfUC = false;
    std::vector<double> energyPrice;        // π_t [$/MWh], Table 2 of the paper
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

// ============================================================
//  CASE 2 – self-UC instance from Morales-Espana et al. (2015)
// ============================================================
ThermalFleet buildFleetPaper()
{
    ThermalFleet f;

    // 24-hour horizon
    f.horizon = 24;

    // 10 generating units  (7 slow-start + 3 quick-start)
    f.unitCount = 10;

    // ------------------------------------------------------------------
    // Unit type flags
    // Units 0-6 : slow-start  (Table 3)
    // Units 7-9 : quick-start (units 8,9,10 from Table 1)
    // ------------------------------------------------------------------
    f.isQuickStart = { false, false, false, false, false, false, false,
                       true,  true,  true };

    // ------------------------------------------------------------------
    // Cost coefficients
    // Table 3 does not publish cost data; we take them from the
    // corresponding quick-start counterparts in Table 1 for units 1-7.
    // Units 8-10 are identical to Table 1 rows 8-10.
    // ------------------------------------------------------------------
    // CLV  [$/MWh]
    f.marginalCost = { 16.19, 17.26, 16.60, 16.50, 19.70, 22.26, 27.74,
                        25.92, 27.74, 27.79 };
    // CNL  [$/h]
    f.fixedCost = { 1000, 970, 700, 680, 450, 370, 480,
                        660, 665, 670 };
    // CSU  [$]
    f.startupCost = { 9000, 10000, 1100, 1120, 1800, 340, 520,
                        60, 60, 60 };
    // CSD  [$]  — paper sets CSD = 0 for all units (Table 1)
    f.shutdownCost = { 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0 };

    // ------------------------------------------------------------------
    // Power limits [MW]  (P, P̲ from Tables 1 & 3)
    // ------------------------------------------------------------------
    f.maximumOutput = { 455, 455, 130, 130, 162, 80, 85,
                        55,  55,  55 };
    f.minimumOutput = { 150, 150,  20,  20,  25, 20, 25,
                        10,  10,  10 };

    // ------------------------------------------------------------------
    // Ramping [MW/h] — not published for the self-UC case; we use
    // generous values (full ramp in 1 h) so ramping is non-binding.
    // ------------------------------------------------------------------
    f.rampUpLimit = { 305, 305, 110, 110, 137, 60, 60,
                        45,  45,  45 };
    f.rampDownLimit = { 305, 305, 110, 110, 137, 60, 60,
                        45,  45,  45 };

    // ------------------------------------------------------------------
    // SU / SD capabilities [MW]  (Table 1 for quick-start; Table 3
    // has SU=SD=Pmin for slow-start units)
    // ------------------------------------------------------------------
    f.startupCap = { 150, 150,  20,  20,  25, 20, 25,
                       25,  25,  25 };
    f.shutdownCap = { 150, 150,  20,  20,  25, 20, 25,
                       25,  25,  25 };

    // ------------------------------------------------------------------
    // Startup / shutdown durations [h]  (SUD, SDD from Table 3)
    // Quick-start units have SUD=SDD=1 (Fig. 3 of the paper)
    // ------------------------------------------------------------------
    f.startupDuration = { 3, 3, 2, 2, 2, 1, 1,  1, 1, 1 };
    f.shutdownDuration = { 2, 2, 2, 2, 2, 1, 1,  1, 1, 1 };

    // ------------------------------------------------------------------
    // Startup/shutdown power trajectories for slow-start units
    // "linear change from 0 to Pmin over SUD hours" (Section 5.1)
    // ------------------------------------------------------------------
    f.startupTraj.resize(10);
    f.shutdownTraj.resize(10);
    for (int g = 0; g < 7; ++g) {
        f.startupTraj[g] = makeStartupTraj(f.minimumOutput[g], f.startupDuration[g]);
        f.shutdownTraj[g] = makeShutdownTraj(f.minimumOutput[g], f.shutdownDuration[g]);
    }
    // quick-start: empty trajectories (handled by the quick-start branch)
    for (int g = 7; g < 10; ++g) {
        f.startupTraj[g] = {};
        f.shutdownTraj[g] = {};
    }

    // ------------------------------------------------------------------
    // Minimum up/down times [h]  (TU = TD from Tables 1 & 3)
    // ------------------------------------------------------------------
    f.minOnline = { 8, 8, 5, 5, 6, 3, 3,  1, 1, 1 };
    f.minOffline = { 8, 8, 5, 5, 6, 3, 3,  1, 1, 1 };

    // ------------------------------------------------------------------
    // Initial conditions (Ste0 and p0 from Table 1; Table 3 mirrors them)
    // initialStatus > 0  means the unit has been ON for that many hours.
    // ------------------------------------------------------------------
    f.initialStatus = { 8, 8, 5, 5, 6, 3, 3,  1, 1, 1 };
    f.initialProduction = { 150, 150, 20, 20, 25, 20, 25,
                            10,  10,  10 };

    // ------------------------------------------------------------------
    // 24-h energy prices [$/MWh]  — Table 2 of the paper
    // ------------------------------------------------------------------
    const std::vector<double> price = {
        13.0,  7.2,  4.6,  3.3,  3.9,  5.9,
         9.8, 15.0, 22.1, 31.3, 33.2, 24.8,
        19.5, 16.3, 14.3, 13.7, 15.0, 17.6,
        20.2, 29.3, 49.5, 53.4, 30.0, 20.2
    };

    // ------------------------------------------------------------------
    // Self-UC mode: price-taker, no external demand (Section 5.1).
    // ------------------------------------------------------------------
    f.selfUC = true;
    f.energyPrice = price;   // π_t from Table 2

    // load and reserve vectors are unused in selfUC mode but kept empty
    f.load.assign(24, 0.0);
    f.reserveUp.assign(24, 0.0);
    f.reserveDown.assign(24, 0.0);

    return f;
}

// ============================================================
//  Helper: solve one UC instance and print results
// ============================================================
void solveInstance(IloEnv& env, const ThermalFleet& fleet, const std::string& caseLabel)
{
    IloModel model(env);

    const int G = fleet.unitCount;
    const int T = fleet.horizon;

    // -------------------------------------------------
    // VARIABLES
    // -------------------------------------------------
    std::vector<std::vector<IloBoolVar>> online(G, std::vector<IloBoolVar>(T));
    std::vector<std::vector<IloBoolVar>> startup(G, std::vector<IloBoolVar>(T));
    std::vector<std::vector<IloBoolVar>> shutdown(G, std::vector<IloBoolVar>(T));

    std::vector<std::vector<IloNumVar>> aboveMinimum(G, std::vector<IloNumVar>(T));
    std::vector<std::vector<IloNumVar>> totalOutput(G, std::vector<IloNumVar>(T));
    std::vector<std::vector<IloNumVar>> producedEnergy(G, std::vector<IloNumVar>(T));
    std::vector<std::vector<IloNumVar>> reservePos(G, std::vector<IloNumVar>(T));
    std::vector<std::vector<IloNumVar>> reserveNeg(G, std::vector<IloNumVar>(T));

    for (int g = 0; g < G; ++g)
        for (int t = 0; t < T; ++t)
        {
            online[g][t] = IloBoolVar(env);
            startup[g][t] = IloBoolVar(env);
            shutdown[g][t] = IloBoolVar(env);

            aboveMinimum[g][t] = IloNumVar(env, 0.0,
                fleet.maximumOutput[g] - fleet.minimumOutput[g], ILOFLOAT);
            totalOutput[g][t] = IloNumVar(env, 0.0, fleet.maximumOutput[g], ILOFLOAT);
            producedEnergy[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            reservePos[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            reserveNeg[g][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
        }

    // -------------------------------------------------
    // OBJECTIVE
    // Cost minimisation (eq. 16/18) for network UC:
    // -------------------------------------------------
    IloExpr objective(env);
    for (int g = 0; g < G; ++g)
    {
        int SUD = fleet.startupDuration[g];
        int SDD = fleet.shutdownDuration[g];
        double CSU_eff = fleet.startupCost[g] + fleet.fixedCost[g] * SUD;
        double CSD_eff = fleet.shutdownCost[g] + fleet.fixedCost[g] * SDD;

        for (int t = 0; t < T; ++t)
        {
            double netMarginal = fleet.marginalCost[g]
                - (fleet.selfUC ? fleet.energyPrice[t] : 0.0);
            objective += netMarginal * producedEnergy[g][t];
            objective += fleet.fixedCost[g] * online[g][t];
            objective += CSU_eff * startup[g][t];
            objective += CSD_eff * shutdown[g][t];
        }
    }
    model.add(IloMinimize(env, objective));

    // -------------------------------------------------
    // (7) LOGIC EQUATIONS
    // -------------------------------------------------
    for (int g = 0; g < G; ++g)
    {
        int u0_prev = (fleet.initialStatus[g] > 0) ? 1 : 0;
        model.add(online[g][0] - u0_prev == startup[g][0] - shutdown[g][0]);
        for (int t = 1; t < T; ++t)
            model.add(online[g][t] - online[g][t - 1] == startup[g][t] - shutdown[g][t]);
    }

    // -------------------------------------------------
    // (8) MINIMUM UP TIME
    // -------------------------------------------------
    for (int g = 0; g < G; ++g)
    {
        int TU = fleet.minOnline[g];
        if (fleet.initialStatus[g] > 0)
        {
            int hoursOn = fleet.initialStatus[g];
            int mustOn = std::max(0, TU - hoursOn);
            for (int t = 0; t < std::min(mustOn, T); ++t)
                model.add(online[g][t] == 1);
        }
        for (int t = TU - 1; t < T; ++t)
        {
            IloExpr lhs(env);
            for (int i = t - TU + 1; i <= t; ++i) lhs += startup[g][i];
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
        if (fleet.initialStatus[g] < 0)
        {
            int hoursOff = -fleet.initialStatus[g];
            int mustOff = std::max(0, TD - hoursOff);
            for (int t = 0; t < std::min(mustOff, T); ++t)
                model.add(online[g][t] == 0);
        }
        for (int t = TD - 1; t < T; ++t)
        {
            IloExpr lhs(env);
            for (int i = t - TD + 1; i <= t; ++i) lhs += shutdown[g][i];
            model.add(lhs <= 1 - online[g][t]);
            lhs.end();
        }
    }

    // -------------------------------------------------
    // (4)(5)(6) GENERATION LIMITS WITH SU/SD
    // -------------------------------------------------
    for (int g = 0; g < G; ++g)
    {
        double Pmax = fleet.maximumOutput[g];
        double Pmin = fleet.minimumOutput[g];
        double SU = fleet.startupCap[g];
        double SD = fleet.shutdownCap[g];

        for (int t = 0; t < T - 1; ++t)
            model.add(
                aboveMinimum[g][t] + reservePos[g][t]
                <= (Pmax - Pmin) * online[g][t]
                - (Pmax - SD) * shutdown[g][t + 1]
                + (SU - Pmin) * startup[g][t + 1]
            );
        model.add(
            aboveMinimum[g][T - 1] + reservePos[g][T - 1]
            <= (Pmax - Pmin) * online[g][T - 1]
        );
    }

    // down-reserve capability
    for (int g = 0; g < G; ++g)
    {
        double Pmin = fleet.minimumOutput[g];
        double SD = fleet.shutdownCap[g];
        for (int t = 0; t < T; ++t)
        {
            model.add(reserveNeg[g][t] <= aboveMinimum[g][t] - (SD - Pmin) * shutdown[g][t]);
            model.add(reserveNeg[g][t] <= aboveMinimum[g][t]);
        }
    }

    // -------------------------------------------------
    // (14)(29) / (12)(30) TOTAL POWER
    // -------------------------------------------------
    for (int g = 0; g < G; ++g)
    {
        double Pmin = fleet.minimumOutput[g];
        for (int t = 0; t < T; ++t)
        {
            IloExpr rhs(env);
            rhs += Pmin * online[g][t];
            if (t + 1 < T) rhs += Pmin * startup[g][t + 1];
            rhs += aboveMinimum[g][t];

            if (!fleet.isQuickStart[g])
            {
                int SUD = fleet.startupDuration[g];
                int SDD = fleet.shutdownDuration[g];
                for (int i = 1; i <= SUD; ++i)
                {
                    int vIdx = t - i + SUD + 2 - 1;
                    double PSU_i = fleet.startupTraj[g][i - 1];
                    if (vIdx >= 0 && vIdx < T) rhs += PSU_i * startup[g][vIdx];
                }
                for (int i = 2; i <= SDD + 1; ++i)
                {
                    int wIdx = t - i + 2 - 1;
                    double PSD_i = fleet.shutdownTraj[g][i - 1];
                    if (wIdx >= 0 && wIdx < T) rhs += PSD_i * shutdown[g][wIdx];
                }
            }
            model.add(totalOutput[g][t] == rhs);
            rhs.end();
        }
    }

    // -------------------------------------------------
    // (15)(31) ENERGY — trapezoidal rule
    // -------------------------------------------------
    for (int g = 0; g < G; ++g)
        for (int t = 0; t < T; ++t)
        {
            IloExpr avg(env);
            avg += totalOutput[g][t];
            avg += (t > 0) ? totalOutput[g][t - 1] : IloExpr(env) + fleet.initialProduction[g];
            model.add(producedEnergy[g][t] == 0.5 * avg);
            avg.end();
        }

    // -------------------------------------------------
    // (27) RAMP-UP
    // -------------------------------------------------
    for (int g = 0; g < G; ++g)
    {
        double Pmin = fleet.minimumOutput[g];
        double RU = fleet.rampUpLimit[g];
        double SU = fleet.startupCap[g];
        int u0_prev = (fleet.initialStatus[g] > 0) ? 1 : 0;
        double p0 = std::max(0.0, fleet.initialProduction[g] - Pmin * u0_prev);

        model.add(aboveMinimum[g][0] + reservePos[g][0] - p0
            <= RU * u0_prev + (SU - Pmin) * startup[g][0]);
        for (int t = 1; t < T; ++t)
            model.add(aboveMinimum[g][t] + reservePos[g][t] - aboveMinimum[g][t - 1]
                <= RU * online[g][t - 1] + (SU - Pmin) * startup[g][t]);
    }

    // -------------------------------------------------
    // (28) RAMP-DOWN
    // -------------------------------------------------
    for (int g = 0; g < G; ++g)
    {
        double Pmin = fleet.minimumOutput[g];
        double RD = fleet.rampDownLimit[g];
        double SD = fleet.shutdownCap[g];
        int u0_prev = (fleet.initialStatus[g] > 0) ? 1 : 0;
        double p0 = std::max(0.0, fleet.initialProduction[g] - Pmin * u0_prev);

        model.add(p0 - (aboveMinimum[g][0] - reserveNeg[g][0])
            <= RD * online[g][0] + (SD - Pmin) * shutdown[g][0]);
        for (int t = 1; t < T; ++t)
            model.add(aboveMinimum[g][t - 1] - (aboveMinimum[g][t] - reserveNeg[g][t])
                <= RD * online[g][t] + (SD - Pmin) * shutdown[g][t]);
    }

    // -------------------------------------------------
    // (19) DEMAND BALANCE — network UC only
    // -------------------------------------------------
    if (!fleet.selfUC)
    {
        for (int t = 0; t < T; ++t)
        {
            IloExpr gen(env);
            for (int g = 0; g < G; ++g) gen += totalOutput[g][t];
            model.add(gen >= fleet.load[t]);
            gen.end();
        }
    }

    // -------------------------------------------------
    // (20)/(21) RESERVES — network UC only
    // -------------------------------------------------
    if (!fleet.selfUC)
    {
        for (int t = 0; t < T; ++t)
        {
            IloExpr res(env);
            for (int g = 0; g < G; ++g) res += reservePos[g][t];
            model.add(res >= fleet.reserveUp[t]);
            res.end();
        }
        for (int t = 0; t < T; ++t)
        {
            IloExpr res(env);
            for (int g = 0; g < G; ++g) res += reserveNeg[g][t];
            model.add(res >= fleet.reserveDown[t]);
            res.end();
        }
    }

    // -------------------------------------------------
    // SOLVE
    // -------------------------------------------------
    IloCplex solver(model);
    solver.setParam(IloCplex::Param::Threads, 1);

    auto sep = [](char c, int n) { return std::string(n, c); };
    auto line = [&]() { std::cout << sep('-', 72) << "\n"; };
    auto dline = [&]() { std::cout << sep('=', 72) << "\n"; };

    std::cout << "\n";
    dline();
    std::cout << "  CASE: " << caseLabel << "\n";
    dline();

    double t0 = solver.getCplexTime();
    if (!solver.solve())
    {
        std::cout << "  No feasible solution found.\n";
        dline();
        objective.end();
        return;
    }
    double t1 = solver.getCplexTime();

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Optimal objective :  " << solver.getObjValue() << " $\n";
    std::cout << "  Best bound        :  " << solver.getBestObjValue() << " $\n";
    std::cout << "  MIP gap           :  " << solver.getMIPRelativeGap() * 100.0 << " %\n";
    std::cout << "  Solve time        :  " << (t1 - t0) << " s\n";
    std::cout << "  B&B nodes         :  " << solver.getNnodes() << "\n";
    std::cout << "  Constraints       :  " << solver.getNrows() << "\n";
    std::cout << "  Variables         :  " << solver.getNcols() << "\n";
    dline();
    std::cout << "\n";

    for (int g = 0; g < G; ++g)
    {
        dline();
        std::cout << "  GENERATOR " << g
            << "  [" << (fleet.isQuickStart[g] ? "quick-start" : "slow-start") << "]"
            << "   P_min=" << fleet.minimumOutput[g] << " MW"
            << "  P_max=" << fleet.maximumOutput[g] << " MW"
            << "  c_var=" << fleet.marginalCost[g] << " $/MWh\n";
        dline();
        std::cout
            << std::left << std::setw(4) << "t"
            << std::setw(14) << "status"
            << std::right << std::setw(10) << "p_abs[MW]"
            << std::setw(10) << "p_tot[MW]"
            << std::setw(9) << "r+[MW]"
            << std::setw(9) << "r-[MW]"
            << std::setw(10) << "E[MWh]"
            << "\n";
        line();
        for (int t = 0; t < T; ++t)
        {
            bool isOn = solver.getValue(online[g][t]) > 0.5;
            bool isSU = solver.getValue(startup[g][t]) > 0.5;
            bool isSD = solver.getValue(shutdown[g][t]) > 0.5;
            std::string st = isOn ? "ON" : "OFF";
            if (isSU) st += "+SU";
            if (isSD) st += "+SD";
            std::cout
                << std::left << std::setw(4) << t
                << std::setw(14) << st
                << std::right << std::setw(10) << solver.getValue(aboveMinimum[g][t])
                << std::setw(10) << solver.getValue(totalOutput[g][t])
                << std::setw(9) << solver.getValue(reservePos[g][t])
                << std::setw(9) << solver.getValue(reserveNeg[g][t])
                << std::setw(10) << solver.getValue(producedEnergy[g][t])
                << "\n";
        }
        int SUD = fleet.startupDuration[g];
        int SDD = fleet.shutdownDuration[g];
        double CSU_eff = fleet.startupCost[g] + fleet.fixedCost[g] * SUD;
        double CSD_eff = fleet.shutdownCost[g] + fleet.fixedCost[g] * SDD;
        double costNL = 0, costSU = 0, costSD = 0, costVar = 0;
        for (int t = 0; t < T; ++t)
        {
            costNL += fleet.fixedCost[g] * solver.getValue(online[g][t]);
            costSU += CSU_eff * solver.getValue(startup[g][t]);
            costSD += CSD_eff * solver.getValue(shutdown[g][t]);
            costVar += fleet.marginalCost[g] * solver.getValue(producedEnergy[g][t]);
        }
        line();
        std::cout << "  Cost:  no-load=" << costNL
            << "  startup=" << costSU
            << "  shutdown=" << costSD
            << "  variable=" << costVar
            << "  total=" << (costNL + costSU + costSD + costVar) << " $\n";
        dline();
        std::cout << "\n";
    }

    // SYSTEM BALANCE CHECK / PROFIT SUMMARY
    dline();
    if (fleet.selfUC)
    {
        std::cout << "  PROFIT SUMMARY (self-UC, price-taker)\n";
        dline();
        std::cout
            << std::left << std::setw(4) << "t"
            << std::right << std::setw(10) << "price[$/MWh]"
            << std::setw(12) << "gen[MW]"
            << std::setw(12) << "energy[MWh]"
            << std::setw(14) << "revenue[$]"
            << "\n";
        line();
        double totalRevenue = 0, totalCostAll = 0;
        for (int t = 0; t < T; ++t)
        {
            double gen = 0, energy = 0;
            for (int g = 0; g < G; ++g)
            {
                gen += solver.getValue(totalOutput[g][t]);
                energy += solver.getValue(producedEnergy[g][t]);
            }
            double rev = fleet.energyPrice[t] * energy;
            totalRevenue += rev;
            std::cout
                << std::left << std::setw(4) << t
                << std::right << std::setw(10) << fleet.energyPrice[t]
                << std::setw(12) << gen
                << std::setw(12) << energy
                << std::setw(14) << rev
                << "\n";
        }
        // total cost = sum over generators
        for (int g = 0; g < G; ++g)
        {
            int SUD = fleet.startupDuration[g];
            int SDD = fleet.shutdownDuration[g];
            double CSU_eff = fleet.startupCost[g] + fleet.fixedCost[g] * SUD;
            double CSD_eff = fleet.shutdownCost[g] + fleet.fixedCost[g] * SDD;
            for (int t = 0; t < T; ++t)
            {
                totalCostAll += fleet.fixedCost[g] * solver.getValue(online[g][t]);
                totalCostAll += CSU_eff * solver.getValue(startup[g][t]);
                totalCostAll += CSD_eff * solver.getValue(shutdown[g][t]);
                totalCostAll += fleet.marginalCost[g] * solver.getValue(producedEnergy[g][t]);
            }
        }
        line();
        std::cout << "  Total revenue : " << totalRevenue << " $\n";
        std::cout << "  Total cost    : " << totalCostAll << " $\n";
        std::cout << "  Net profit    : " << (totalRevenue - totalCostAll) << " $\n";
        dline();
    }
    else
    {
        std::cout
            << std::left << std::setw(4) << "t"
            << std::right << std::setw(10) << "demand[MW]"
            << std::setw(11) << "gen[MW]"
            << std::setw(8) << "slack"
            << std::setw(9) << "r+req"
            << std::setw(9) << "sum_r+"
            << std::setw(9) << "r-req"
            << std::setw(9) << "sum_r-"
            << "\n";
        line();
        bool allOk = true;
        for (int t = 0; t < T; ++t)
        {
            double gen = 0, sumRp = 0, sumRn = 0;
            for (int g = 0; g < G; ++g)
            {
                gen += solver.getValue(totalOutput[g][t]);
                sumRp += solver.getValue(reservePos[g][t]);
                sumRn += solver.getValue(reserveNeg[g][t]);
            }
            double slack = gen - fleet.load[t];
            bool ok = (slack >= -1e-4)
                && (sumRp >= fleet.reserveUp[t] - 1e-4)
                && (sumRn >= fleet.reserveDown[t] - 1e-4);
            if (!ok) allOk = false;
            std::cout
                << std::left << std::setw(4) << t
                << std::right << std::setw(10) << fleet.load[t]
                << std::setw(11) << gen
                << std::setw(8) << (ok ? "OK" : "FAIL")
                << std::setw(9) << fleet.reserveUp[t]
                << std::setw(9) << sumRp
                << std::setw(9) << fleet.reserveDown[t]
                << std::setw(9) << sumRn
                << "\n";
        }
        line();
        std::cout << "  Overall: "
            << (allOk ? "ALL CONSTRAINTS SATISFIED"
                : "WARNING — constraint violation detected")
            << "\n";
        dline();
    }

    objective.end();
}

int main()
{
    IloEnv env;

    try
    {
        // Case 1: original custom instance
        solveInstance(env, buildFleet(), "CASE 1 — Custom test instance (T=6h, 3 units)");

        // Case 2: paper instance
        solveInstance(env, buildFleetPaper(), "CASE 2 — Morales-Espana et al. (2015), Table 1+3, T=24h, 10 units");
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