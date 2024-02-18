#include <cstring>
#include <ilcplex/ilocplex.h>
#include <string>

typedef IloArray<IloNumVarArray> NumVarMatrix;

using namespace std;

int main(int argc, char *argv[]){
	IloEnv env;

	try{
		const IloInt n = 9, m = 2, d = 5, T = 200, H = 10;
		IloNumArray s(env, m, 2, 3);
		IloNumArray r(env, 1, 4);
		IloNumArray L(env, n, 18, 14, 12, 6, 36, 9, 9, 12, 12);
		IloNumArray C(env, n,  4,  4,  6, 6, 12, 6, 2,  2,  4);
		IloNumArray h(env, n,  3,  3,  4, 4,  8, 6, 2,  2,  4);

		IloModel model(env);
		IloNumVar z(env, 0., IloInfinity, ILOFLOAT, "z");
		IloArray<NumVarMatrix> x(env, m);
		IloArray<NumVarMatrix> f(env, m);
		IloArray<NumVarMatrix> p(env, m);
		char name[8]; strncpy(name, "xij\0", 8);
		char namef[8]; strncpy(namef, "fij\0", 8);
		char namep[8]; strncpy(namep, "pijl\0", 8);
		for (int i = 0; i < m; i++) {
			x[i] = NumVarMatrix(env, n);
			f[i] = NumVarMatrix(env, n);
			p[i] = NumVarMatrix(env, n+1);
			for (int j = 0; j < n; j++){
				x[i][j] = IloNumVarArray(env, T-L[j]-s[i]+r[0], 0., 1., ILOBOOL);
				f[i][j] = IloNumVarArray(env, T, 0., 1., ILOBOOL);
				name[1] = i + '0';
				name[2] = j + '0';
				x[i][j].setNames(name);
				namef[1] = i + '0';
				namef[2] = j + '0';
				f[i][j].setNames(namef);
			}
			for (int j = 0; j < n+1; j++){
				p[i][j] = IloNumVarArray(env, n+1, 0., 1., ILOBOOL);
				namep[1] = i + '0';
				namep[2] = j + '0';
				for (int l = 0; l < n+1; l++) {
					namep[3] = l + '0';
					p[i][j][l].setName(namep);
				}
				model.add(p[i][j][j] == 0);
			}
		}

		model.add(IloMinimize(env, z));
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++){
				IloExpr exp(env);
				for (int t = 0; t < T; t++){
					exp += (t + r[0])*f[i][j][t];
				}
				exp +=  C[j];
				model.add(z >= exp);
			}
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++){
				for (int t = 0; t < T-L[j]-s[i]+r[0]; t++){
					IloExpr exp(env);
					for (int k = t+L[j]+s[i]-r[0]; k < T; k++) {
						exp += f[i][j][k];
					}
					model.add(x[i][j][t] <= exp);
				}
			}
		}

		for (int j = 0; j < n; j++){
			IloExpr exp(env);
			for (int i = 0; i < m; i++) {
				for (int t = 0; t < T-L[j]-s[i]+r[0]; t++){
					exp += x[i][j][t];
				}
			}
			model.add(exp >= 1);
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++){
				IloExpr exp(env);
				for (int t = 0; t < T; t++){
					exp += f[i][j][t];
				}
				model.add(exp <= 1);
			}
		}

		for (int i = 0; i < m; i++) {
			for (int t = 0; t < T; t++) {
				IloExpr exp(env);
				for (int j = 0; j < n; j++) {
					for (int k = 0; k <= t; k++) {
						if (t-L[j]-s[i] <= k && k < T-L[j]-s[i]+r[0]) {
							exp += x[i][j][k];
						}
					}
				}
				model.add(exp <= 1);
			}
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int l = 0; l < n; l++) {
					if (l != i){
						IloExpr exp(env);
						for (int k = 0; k < T-L[j]-s[i]+r[0]; k++) {
							exp += (k + s[i])*x[i][j][k];
						}
						for (int k = 0; k < T-L[l]-s[i]+r[0]; k++) {
							exp += (-k)*x[i][l][k];
						}
						exp += (-T - s[0] -L[j])*(1 - p[i][j][l]) + L[j];
						model.add(exp <= 0);
					}
				}
			}
		}

		for (int j = 0; j < n; j++){
			IloExpr exp(env);
			for (int i = 0; i < m; i++) {
				for (int t = 0; t < T; t++){
					exp += f[i][j][t];
				}
			}
			model.add(exp <= 1);
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int l = 0; l < n; l++) {
					if (l != i){
					IloExpr exp(env);
						for (int k = 0; k < T; k++) {
							exp += (k + r[0])*f[i][j][k];
							exp += (-k)*f[i][l][k];
						}
						exp += (-T - r[0] - C[j])*(1 - p[i][j][l]) + C[j];
						model.add(exp <= 0);
					}
				}
			}
		}

		for (int u = 0; u < m; u++) {
			for (int v = u+1; v < m; v++) {
				for (int t = 0; t < T; t++) {
					IloExpr exp(env);
					for (int j = 0; j < n; j++) {
						for (int k = t-C[j]-r[0]-d; k <= t; k++) {
							if (k >= 0) {
								exp += f[u][j][k];
								exp += f[v][j][k];
							}
						}
					}
					model.add(exp <= 1);
				}
			}
		}

		for (int i = 0; i < m; i++) {
			for (int t = 0; t < T; t++) {
				IloExpr exp(env);
				for (int j = 0; j < n; j++) {
					for (int k = 0; k <= t-s[i]; k++) {
						if (k < T-L[j]-s[i]+r[0]) exp += h[j]*x[i][j][k];
					}
					for (int k = 0; k <= t-r[0]-C[j]; k++) {
						if (k < T-L[j]-s[i]+r[0]) exp += -h[j]*f[i][j][k];
					}
				}
				model.add(exp <= H);
			}
		}


		///
		// for (int i = 0; i < m; i++) {
		// 	for (int t = 0; t < T; t++) {
		// 		IloExpr exp(env);
		// 		for (int j = 0; j < n; j++) {
		// 			for (int k = 0; k <= t; k++) {
		// 				if (t-L[j]-s[i] <= k && k < T-L[j]-s[i]+r[0]) {
		// 					exp += x[i][j][k];
		// 				}
		// 			}
		// 		}
		// 		model.add(exp <= 1);
		// 	}
		// }
		///

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				IloExpr exp(env);
				for (int t = 0; t < T-L[j]-s[i]+r[0]; t++) {
					exp += x[i][j][t];
				}
				for (int l = 0; l < n+1; l++) {
					if (l != j){
						exp += (-1.)*p[i][j][l];
					}
				}
				model.add(exp == 0);
			}
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				IloExpr exp(env);
				for (int t = 0; t < T-L[j]-s[i]+r[0]; t++) {
					exp += x[i][j][t];
				}
				for (int l = 0; l < n+1; l++) {
					if (l != j){
						exp += (-1.)*p[i][l][j];
					}
				}
				model.add(exp == 0);
			}
		}

		for (int i = 0; i < m; i++) {
			IloExpr exp(env);
			for (int j = 0; j < n; j++) {
				exp += p[i][n][j];
			}
			model.add(exp == 1);
		}

		for (int i = 0; i < m; i++) {
			IloExpr exp(env);
			for (int j = 0; j < n; j++) {
				exp += p[i][j][n];
			}
			model.add(exp == 1);
		}
				
		IloCplex cplex(model);
		cplex.exportModel("out.lp");

		cplex.setParam(IloCplex::Param::DetTimeLimit, 360000);
		int err = cplex.solve();
		if (!err) throw IloAlgorithm::Exception("Failed to optimize.");

		cplex.out() << endl;
		cplex.out() << "Solution: " << cplex.getStatus() << '\n';
		cplex.out() << "Objective value = " << cplex.getObjValue() << '\n';
		cplex.out() << endl;

		for (int i = 0; i < m; i++) {
			cplex.out() << "Machine " << i+1 << ":...";
			for (int t = 0; t < T; t++){
				bool idle = true;
				for (int j = 0; j < n; j++){
					if (t < T-L[j]-s[i]+r[0] && cplex.getValue(x[i][j][t]) > 0.1) {
						idle = false;
						for (int k = 0; k < s[i]; k++) {
							cplex.out() << 's';
						}
						for (int k = 0; k < L[j]; k++) {
							cplex.out() <<  to_string(j);
						}
						t += L[j] + s[i];
						break;
					}
				}
				if (idle) cplex.out() << 'x';
			}
			cplex.out() << endl << "Cutting " << i+1 << ":...";
			for (int t = 0; t < T; t++){
				bool idle = true;
				for (int j = 0; j < n; j++){
					if (cplex.getValue(f[i][j][t]) > 0.1) {
						idle = false;
						for (int k = 0; k < r[0]; k++) {
							cplex.out() << 's';
						}
						for (int k = 0; k < C[j]; k++) {
							cplex.out() <<  to_string(j);
						}
						t += C[j] + r[0];
						break;
					}
				}
				if (idle) cplex.out() << 'x';
			}
			cplex.out() << "\n\n";
		}
		cplex.out() << "\nx:  ";
		for (int t = 0; t < T; t++){
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++){
					if (t < T-L[j]-s[i]+r[0] && cplex.getValue(x[i][j][t]) > 0.1) {
						cplex.out() << "  x[" << i << "][" <<  to_string(j) << "][" << to_string(t) << "]";
					}
				}
			}
		}
		cplex.out() << "\nf:  ";
		for (int t = 0; t < T; t++){
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++){
					if (cplex.getValue(f[i][j][t]) > 0.1) {
						cplex.out() << "  f[" << i << "][" <<  to_string(j) << "][" << to_string(t) << "]";
					}
				}
			}
		}
		cplex.out() << "\np:  ";
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n+1; j++){
				for (int l = 0; l < n+1; l++){
					if (cplex.getValue(p[i][j][l]) > 0.1) {
						cplex.out() << "  p[" << i << "][" <<  to_string(j) << "][" << to_string(l) << "]";
					}
				}
			}
		}
		cplex.out() << endl;

		// cplex.out() << endl;
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
