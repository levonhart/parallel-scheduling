#include <cstring>
#include <ilcplex/ilocplex.h>
#include <string>

typedef IloArray<IloNumVarArray> NumVarMatrix;

using namespace std;

int main(int argc, char *argv[]){
	IloEnv env;

	try {
		const IloInt n = 4, m = 2;
		IloNumArray s(env, m, 2, 3);
		IloNumArray L(env, n, 12, 10, 8, 4);

		IloModel model(env);
		IloNumVar z(env, 0., IloInfinity, ILOFLOAT, "z");
		NumVarMatrix y(env, m);
		char name[8]; strncpy(name, "y(i)\0", 8);
		for (int i = 0; i < m; i++) {
			y[i] = IloNumVarArray(env, n, 0., 1., ILOBOOL);
			name[2] = i + '0';
			y[i].setNames(name);
		}

		model.add(IloMinimize(env, z));
		for (int i = 0; i < m; i++) {
		IloExpr exp(env);
			for (int j = 0; j < n; j++) {
				exp += (s[i] + L[j])*y[i][j];
			}
			model.add(z >= exp);
		}

		for (int j = 0; j < n; j++) {
			IloExpr exp(env);
			for (int i = 0; i < m; i++) {
				exp += y[i][j];
			}
			model.add(exp >= 1);
		}
				
		IloCplex cplex(model);
		cplex.exportModel("out.lp");

		cplex.setParam(IloCplex::Param::DetTimeLimit, 3600);
		int err = cplex.solve();
		if (!err) throw IloAlgorithm::Exception("Failed to optimize");

		cplex.out() << endl;
		cplex.out() << "Solution: " << cplex.getStatus() << '\n';
		cplex.out() << "Objective value = " << cplex.getObjValue() << '\n';
		cplex.out() << endl;

		for (int i = 0; i < m; i++) {
			cplex.out() << "Machine " << i << " jobs: ";
			for (int j = 0; j < n; j++) {
				cplex.out() << " " << (cplex.getValue(y[i][j]) > 0.1 ? to_string(j) : " ");
			}
			cplex.out() << '\n';
		}
		cplex.out() << endl;

	} catch (IloException& ex) {
		cerr << "Concert exception caught" << endl;
		env.end();
		return EXIT_FAILURE;
	} catch (...) {
		cerr << "Unknown exception caught" << endl;
		env.end();
		return EXIT_FAILURE;
	}

	env.end();

	return 0;
}
