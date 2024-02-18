#include <cstring>
#include <ilcplex/ilocplex.h>
#include <string>

typedef IloArray<IloNumVarArray> NumVarMatrix;

using namespace std;

int main(int argc, char *argv[]){
	IloEnv env;

	try{
		const IloInt n = 4, m = 1, T = 114;
		IloNumArray s(env, m, 2);
		IloNumArray r(env, m, 4);
		IloNumArray L(env, n, 18, 14, 12, 6);
		IloNumArray C(env, n, 4, 4, 6, 6);

		IloModel model(env);
		IloNumVar z(env, 0., IloInfinity, ILOFLOAT, "z");
		NumVarMatrix x(env, n);
		NumVarMatrix f(env, n);
		char name[8]; strncpy(name, "x1i\0", 8);
		char namef[8]; strncpy(namef, "f1i\0", 8);
		for (int j = 0; j < n; j++){
			x[j] = IloNumVarArray(env, T-L[j]-s[0]+r[0], 0., 1., ILOBOOL);
			f[j] = IloNumVarArray(env, T, 0., 1., ILOBOOL);
			name[2] = j + '0';
			x[j].setNames(name);
			namef[2] = j + '0';
			f[j].setNames(namef);
		}

		model.add(IloMinimize(env, z));
		for (int j = 0; j < n; j++){
			IloExpr exp(env);
			for (int t = 0; t < T; t++){
				exp += (t + r[0])*f[j][t];
			}
			exp +=  C[j];
			model.add(z >= exp);
		}

		for (int j = 0; j < n; j++){
			for (int t = 0; t < T-L[j]-s[0]+r[0]; t++){
				IloExpr exp(env);
				for (int k = t+L[j]+s[0]-r[0]; k < T; k++) {
					exp += f[j][k];
				}
				model.add(x[j][t] <= exp);
			}
		}

		for (int j = 0; j < n; j++){
			IloExpr exp(env);
			for (int t = 0; t < T-L[j]-s[0]+r[0]; t++){
				exp += x[j][t];
			}
			model.add(exp == 1);
		}

		for (int t = 0; t < T; t++) {
			IloExpr exp(env);
			for (int j = 0; j < n; j++) {
				for (int k = 0; k <= t; k++) {
					if (t-L[j]-s[0] <= k && k < T-L[j]-s[0]+r[0]) {
						exp += x[j][k];
					}
				}
			}
			model.add(exp <= 1);
		}

		for (int j = 0; j < n; j++){
			IloExpr exp(env);
			for (int t = 0; t < T; t++){
				exp += f[j][t];
			}
			model.add(exp <= 1);
		}

		for (int t = 0; t < T; t++) {
			for (int j = 0; j < n; j++) {
				for (int l = j + 1; l < n; l++) {
					IloExpr exp(env);
					for (int k = 0; k <= t; k++) {
						if (t-C[j]-r[0] <= k) {
							exp += f[j][k];
						}
						if (t-C[l]-r[0] <= k) {
							exp += f[l][k];
						}
					}

					model.add(exp <= 1);
				}
			}
		}
				
		IloCplex cplex(model);
		// cplex.exportModel("out.lp");

		cplex.setParam(IloCplex::Param::DetTimeLimit, 3600);
		int err = cplex.solve();
		if (!err) throw IloAlgorithm::Exception("Failed to optimize.");

		cplex.out() << endl;
		cplex.out() << "Solution: " << cplex.getStatus() << '\n';
		cplex.out() << "Objective value = " << cplex.getObjValue() << '\n';
		cplex.out() << endl;

		cplex.out() << "Machine 1:...";
		for (int t = 0; t < T; t++){
			bool idle = true;
			for (int j = 0; j < n; j++){
				if (t < T-L[j]-s[0]+r[0] && cplex.getValue(x[j][t]) > 0.1) {
					idle = false;
					for (int i = 0; i < s[0]; i++) {
						cplex.out() << 's';
					}
					for (int i = 0; i < L[j]; i++) {
						cplex.out() <<  to_string(j);
					}
					t += L[j] + s[0];
					break;
				}
			}
			if (idle) cplex.out() << 'x';
		}
		cplex.out() << endl << "Cutting:.....";
		for (int t = 0; t < T; t++){
			bool idle = true;
			for (int j = 0; j < n; j++){
				if (cplex.getValue(f[j][t]) > 0.1) {
					idle = false;
					for (int i = 0; i < r[0]; i++) {
						cplex.out() << 's';
					}
					for (int i = 0; i < C[j]; i++) {
						cplex.out() <<  to_string(j);
					}
					t += C[j] + r[0];
					break;
				}
			}
			if (idle) cplex.out() << 'x';
		}
		cplex.out() << "\nx:  ";
		for (int t = 0; t < T; t++){
			for (int j = 0; j < n; j++){
				if (t < T-L[j]-s[0]+r[0] && cplex.getValue(x[j][t]) > 0.1) {
					cplex.out() << "  x1[" <<  to_string(j) << "][" << to_string(t) << "]";
				}
			}
		}
		cplex.out() << "\nf:  ";
		for (int t = 0; t < T; t++){
			for (int j = 0; j < n; j++){
				if (cplex.getValue(f[j][t]) > 0.1) {
					cplex.out() << "  f1[" <<  to_string(j) << "][" << to_string(t) << "]";
				}
			}
		}
		cplex.out() << endl;
		// IloNumArray vals(env);
		// for (int j = 0; j < n; j++) {
		// 	cplex.getValues(vals, x[j]);
		// 	env.out() << "x[" << to_string(j) << "] =\n" << vals << "\n\n";
		// }
		// for (int j = 0; j < n; j++) {
		// 	cplex.getValues(vals, f[j]);
		// 	env.out() << "f[" << to_string(j) << "] =\n" << vals << "\n\n";
		// }

	} catch (IloException& ex){
		cerr << "Concert exception: " << ex << endl;
		env.end();
		return EXIT_FAILURE;
	} catch (...){
		cerr << "Unknown exception." << endl;
		env.end();
		return EXIT_FAILURE;
	}

	env.end();

	return 0;
}
