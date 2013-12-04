/*
1D ADM
2D: Burger at least
2D ADM
2D unstructured
*/

var panel;
var canvas;
var waveVtxBuf, waveStateBuf;
var xmin = 0;
var xmax = 100; 
var ymin = -10;
var ymax = 50;
var gridstep = 10;
var useCFL = true;
var fixedDT = .2;
var gaussSeidelIterations = 20;

function isnan(x) {
	return x != x;
}

function mat33invert(out, a) {
	var det = a[0][0] * a[1][1] * a[2][2]
			+ a[1][0] * a[2][1] * a[0][2]
			+ a[2][0] * a[0][1] * a[1][2]
			- a[2][0] * a[1][1] * a[0][2]
			- a[1][0] * a[0][1] * a[2][2]
			- a[0][0] * a[2][1] * a[1][2];
	if (det == 0) {
		for (var j = 0; j < 3; ++j) {
			for (var i = 0; i < 3; ++i) {
				console.log('a('+i+','+j+') = '+a[j][i]);
			}
		}
		throw 'singular!';
	}
	var invDet = 1 / det;
	for (var j = 0; j < 3; ++j) {
		var j1 = (j + 1) % 3;
		var j2 = (j + 2) % 3;
		for (var i = 0; i < 3; ++i) {
			var i1 = (i + 1) % 3;
			var i2 = (i + 2) % 3;
			out[i][j] = invDet * (a[j1][i1] * a[j2][i2] - a[j1][i2] * a[j2][i1]);
		}
	}
}

/*
output:
matrix
eigenvalues
eigenvectors

input:
velocity
hTotal
speedOfSound
*/
function buildEigenstate(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocity, hTotal, gamma) {
	//calculate matrix & eigenvalues & vectors at interface from state at interface
	var speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * velocity * velocity));	
	//matrix, listed per column
	matrix[0][0] = 0;
	matrix[0][1] = (gamma - 3) / 2 * velocity * velocity;
	matrix[0][2] = velocity * ((gamma - 1) / 2 * velocity * velocity - hTotal);
	matrix[1][0] = 1;
	matrix[1][1] = (3 - gamma) * velocity;
	matrix[1][2] = hTotal - (gamma - 1) * velocity * velocity;
	matrix[2][0] = 0;
	matrix[2][1] = gamma - 1;
	matrix[2][2] = gamma * velocity;

	//eigenvalues: min, mid, max
	eigenvalues[0] = velocity - speedOfSound;
	eigenvalues[1] = velocity;
	eigenvalues[2] = velocity + speedOfSound;
	//min eigenvector
	eigenvectors[0][0] = 1;
	eigenvectors[0][1] = velocity - speedOfSound;
	eigenvectors[0][2] = hTotal - speedOfSound * velocity;
	//mid eigenvector
	eigenvectors[1][0] = 1;
	eigenvectors[1][1] = velocity;
	eigenvectors[1][2] = .5 * velocity * velocity;
	//max eigenvector
	eigenvectors[2][0] = 1;
	eigenvectors[2][1] = velocity + speedOfSound;
	eigenvectors[2][2] = hTotal + speedOfSound * velocity;
	//calculate eigenvector inverses ... 
	mat33invert(eigenvectorsInverse, eigenvectors);
}

var fluxMethods = {
	'donor cell' : function(r) { return 0; },
	'Lax-Wendroff' : function(r) { return 1; },
	
	//these two are no good with the Godunov (Riemann) solver
	'Beam-Warming' : function(r) { return r; },
	'Fromm' : function(r) { return .5 * (1 + r); },

	//Wikipedia
	CHARM : function(r) { if (r < 0) return 0; return r*(3*r+1)/((r+1)*(r+1)); },
	HCUS : function(r) { return Math.max(0, 1.5 * (r + Math.abs(r)) / (r + 2) ); },
	HQUICK : function(r) { return Math.max(0, 2 * (r + Math.abs(r)) / (r + 3) ); },
	Koren : function(r) { return Math.max(0, Math.min(2*r, (1 + 2*r)/3 ,2) ); },
	minmod : function(r) { return Math.max(0, Math.min(r,1) ); },
	Oshker : function(r) { return Math.max(0, Math.min(r,1.5) ); },	//replace 1.5 with 1 <= beta <= 2	
	ospre : function(r) { return .5 * (r*r + r) / (r*r + r + 1); },
	smart : function(r) { return Math.max(0, Math.min(2 * r, .25 + .75 * r, 4)); },
	Sweby : function(r) { return Math.max(0, Math.min(1.5 * r, 1), Math.min(r, 1.5)); },	//replace 1.5 with 1 <= beta <= 2
	UMIST : function(r) { return Math.max(0, Math.min(2*r, .75 + .25*r, .25 + .75*r, 2)); },	
	'van Albada 1' : function(r) { return (r * r + r) / (r * r + 1); },
	'van Albada 2' : function(r) { return 2 * r / (r * r + 1); },
	
	'van Leer' : function(r) { return (r + Math.abs(r)) / (1 + Math.abs(r)); },
	'monotonized central' : function(r) { return Math.max(0, Math.min(2, .5 * (1 + r), 2 * r)); },
	superbee : function(r) { return Math.max(0,Math.min(1,2*r),Math.min(2,r)); }
};

var boundaryMethods = {
	periodic : function(nx,q) {
		q[0][0] = q[nx-4][0];
		q[1][0] = q[nx-3][0];
		q[nx-2][0] = q[2][0];
		q[nx-1][0] = q[3][0];
		q[0][1] = q[nx-4][1];
		q[1][1] = q[nx-3][1];
		q[nx-2][1] = q[2][1];
		q[nx-1][1] = q[3][1];
		q[0][2] = q[nx-4][2];
		q[1][2] = q[nx-3][2];
		q[nx-2][2] = q[2][2];
		q[nx-1][2] = q[3][2];
	},
	mirror : function(nx,q) {
		q[0][0] = q[3][0];
		q[1][0] = q[2][0];
		q[nx-2][0] = q[nx-3][0];
		q[nx-1][0] = q[nx-4][0];
		q[0][1] = -q[3][1];
		q[1][1] = -q[2][1];
		q[nx-2][1] = -q[nx-3][1];
		q[nx-1][1] = -q[nx-4][1];
		q[0][2] = q[3][2];
		q[1][2] = q[2][2];
		q[nx-2][2] = q[nx-3][2];
		q[nx-1][2] = q[nx-4][2];
	},
	dirichlet : function(nx,q) {
		q[0][0] = 0;
		q[1][0] = 0;
		q[nx-2][0] = 0;
		q[nx-1][0] = 0;
		q[0][1] = 0;
		q[1][1] = 0;
		q[nx-2][1] = 0;
		q[nx-1][1] = 0;
		q[0][2] = 0;
		q[1][2] = 0;
		q[nx-2][2] = 0;
		q[nx-1][2] = 0;
	},
	constant : function(nx,q) {
		q[0][0] = q[1][0] = q[2][0];
		q[nx-1][0] = q[nx-2][0] = q[nx-3][0];
		q[0][1] = q[1][1] = q[2][1];
		q[nx-1][1] = q[nx-2][1] = q[nx-3][1];
		q[0][2] = q[1][2] = q[2][2];
		q[nx-1][2] = q[nx-2][2] = q[nx-3][2];
	}
};

//'this' is HydroState
var integrationMethods = {
	'Forward Euler' : {
		initPressure : function() {
			for (var i = 0; i < this.nx; ++i) {
				var u = this.q[i][1] / this.q[i][0];
				var energyTotal = this.q[i][2] / this.q[i][0];
				var energyKinematic = .5 * u * u;
				var energyThermal = energyTotal - energyKinematic;
				this.pressure[i] = (this.gamma - 1) * this.q[i][0] * energyThermal;
			}
		},
		applyMomentumDiffusion : function(dt) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				this.q[i][1] -= dt * (this.pressure[i+1] - this.pressure[i-1]) / (this.x[i+1] - this.x[i-1]);
			}
		},
		applyWorkDiffusion : function(dt) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				var u_inext = this.q[i+1][1] / this.q[i+1][0];
				var u_iprev = this.q[i-1][1] / this.q[i-1][0];
				this.q[i][2] -= dt * (this.pressure[i+1] * u_inext - this.pressure[i-1] * u_iprev) / (this.x[i+1] - this.x[i-1]);
			}
		},
		advect : {
			Burgers : function(dt) {
				assert(this.x.length == this.nx);
				assert(this.xi.length == this.nx + 1);
				assert(this.q.length == this.nx);
				assert(this.ui.length == this.nx + 1);
			
				//get velocity at interfaces from state
				for (var ix = this.nghost-1; ix < this.nx+this.nghost-2; ++ix) {
					this.ui[ix] = .5 * (this.q[ix][1] / this.q[ix][0] + this.q[ix-1][1] / this.q[ix-1][0]);
				}
				this.ui[0] = this.ui[this.nx] = 0;

				//compute flux and advect for each state vector
				for (var j = 0; j < 3; ++j) {
					//r_{i-1/2} flux limiter
					for (var i = this.nghost; i < this.nx+this.nghost-3; ++i) {
						var dq = this.q[i][j] - this.q[i-1][j];
						if (Math.abs(dq) > 0) {
							if (this.ui[i] >= 0) {
								this.r[i][j] = (this.q[i-1][j] - this.q[i-2][j]) / dq;
							} else {
								this.r[i][j] = (this.q[i+1][j] - this.q[i][j]) / dq;
							}
						} else {
							this.r[i][j] = 0;
						}
					}
					this.r[0][j] = this.r[1][j] = this.r[this.nx-1][j] = this.r[this.nx][j] = 0;

					//construct flux:
					for (var i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
						//flux limiter
						var phi = fluxMethods[this.fluxMethod](this.r[i][j]);
						if (this.ui[i] >= 0) {
							this.flux[i][j] = this.ui[i] * this.q[i-1][j];
						} else {
							this.flux[i][j] = this.ui[i] * this.q[i][j];
						}
						var delta = phi * (this.q[i][j] - this.q[i-1][j]);
						var dx = this.x[i] - this.x[i-1];
						this.flux[i][j] += delta * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
					}
					this.flux[0][j] = this.flux[this.nx][j] = 0;

					//update cells
					for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
						this.q[i][j] -= dt * (this.flux[i+1][j] - this.flux[i][j]) / (this.xi[i+1] - this.xi[i]);
					}
				}			
			},
			
			/*
			relation:
			
			eigenvalues:
			lambda 1 = u - Cs
			lambda 2 = u
			lambda 3 = u + Cs
			eigenvectors:
			e 1 = (1, u - Cs, h_total - Cs u)
			e 2 = (1, u, .5*u^2)
			e 3 = (1, u + Cs, h_total + Cs u)
			rho = q0
			u = q1 / q0
			h_total = e_total + P / rho
			Cs = sqrt(gamma P / rho)
			*/
			'Riemann / Roe' : function(dt) {
				for (var ix = 1; ix < this.nx; ++ix) {
					for (var j = 0; j < 3; ++j) {
						this.interfaceDeltaQTilde[ix][j] = this.interfaceEigenvectorsInverse[ix][0][j] * (this.q[ix][0] - this.q[ix-1][0])
												+ this.interfaceEigenvectorsInverse[ix][1][j] * (this.q[ix][1] - this.q[ix-1][1])
												+ this.interfaceEigenvectorsInverse[ix][2][j] * (this.q[ix][2] - this.q[ix-1][2]);
					}
				}
				
				for (var j = 0; j < 3; ++j) {
					this.interfaceDeltaQTilde[0][j] = 0;
					this.interfaceDeltaQTilde[this.nx][j] = 0;
				}
				
				for (var ix = this.nghost; ix < this.nx+this.nghost-3; ++ix) {
					for (var j = 0; j < 3; ++j) {
						var interfaceDeltaQTilde = this.interfaceDeltaQTilde[ix][j];
						if (Math.abs(interfaceDeltaQTilde) > 0) {
							if (this.interfaceEigenvalues[j] > 0) {
								this.rTilde[ix][j] = this.interfaceDeltaQTilde[ix-1][j] / interfaceDeltaQTilde;
							} else {
								this.rTilde[ix][j] = this.interfaceDeltaQTilde[ix+1][j] / interfaceDeltaQTilde;
							}
						} else {
							this.rTilde[ix][j] = 0;
						}
					}
				}

				//..and keep the boundary r's zero	
				for (var j = 0; j < 3; ++j) {
					this.rTilde[0][j] = this.rTilde[1][j] = this.rTilde[this.nx-1][j] = this.rTilde[this.nx][j] = 0;
				}
				/*
				for (var ix = 0; ix < this.nghost; ++ix) {
					for (var j = 0; j < 3; ++j) {
						this.rTilde[ix][j] = 0;
						this.rTilde[this.nx-ix][j] = 0;
					}	
				}
				*/
				
				//transform cell q's into cell qTilde's (eigenspace)
				// ... so q_{i-1/2}L = q_{i-1}, q_{i-1/2}R = q_i
				// qTilde_{i-1/2}L = E_{i-1/2}^-1 q_{i-1}, qTilde_{i-1/2}R = E_{i-1/2}^-1 q_i
				//use them to detemine qTilde's at boundaries
				//use them (and eigenvalues at boundaries) to determine fTilde's at boundaries
				//use them (and eigenvectors at boundaries) to determine f's at boundaries
				//use them to advect, like good old fluxes advect
				var fluxTilde = [];
				var fluxAvg = [];

				//qi[ix] = q_{i-1/2} lies between q_{i-1} = q[i-1] and q_i = q[i]
				//(i.e. qi[ix] is between q[ix-1] and q[ix])
				//Looks good according to "Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro, p.191
				for (var ix = 1; ix < this.nx; ++ix) {
					//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
					//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
					for (var j = 0; j < 3; ++j) {
						fluxAvg[j] = .5 * ( 
							this.interfaceMatrix[ix][0][j] * (this.q[ix-1][0] + this.q[ix][0])
							+ this.interfaceMatrix[ix][1][j] * (this.q[ix-1][1] + this.q[ix][1])
							+ this.interfaceMatrix[ix][2][j] * (this.q[ix-1][2] + this.q[ix][2]));
					}

					//calculate flux
					for (var j = 0; j < 3; ++j) {
						var theta = 0;
						if (this.interfaceEigenvalues[ix][j] >= 0) {
							theta = 1;
						} else {
							theta = -1;
						}
						
						var phi = fluxMethods[this.fluxMethod](this.rTilde[ix][j]);
						var dx = this.xi[ix] - this.xi[ix-1];
						var epsilon = this.interfaceEigenvalues[ix][j] * dt / dx;
						
						//flux[ix][k] = fluxTilde[ix][j] * interfaceEigenvectors[ix][k][j]
						//flux in eigenvector basis is the q vector transformed by the inverse then scaled by the eigenvalue
						//should the eigenvalue be incorperated here, after flux limiter is taken into account, or beforehand?
						//1D says after, but notes say before ...
						var deltaFluxTilde = this.interfaceEigenvalues[ix][j] * this.interfaceDeltaQTilde[ix][j];
						
						fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
					}

					//reproject fluxTilde back into q
					for (var j = 0; j < 3; ++j) {
						this.flux[ix][j] = fluxAvg[j]
							+ this.interfaceEigenvectors[ix][0][j] * fluxTilde[0] 
							+ this.interfaceEigenvectors[ix][1][j] * fluxTilde[1] 
							+ this.interfaceEigenvectors[ix][2][j] * fluxTilde[2];
					}

				}
				
				//zero boundary flux
				for (var j = 0; j < 3; ++j) {
					this.flux[0][j] = this.flux[this.nx][j] = 0;
				}

				//update cells
				for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
					for (var j = 0; j < 3; ++j) {
						this.q[i][j] -= dt * (this.flux[i+1][j] - this.flux[i][j]) / (this.xi[i+1] - this.xi[i]);
					}
				}
			}
		}
	},
	//TODO 3x3 block tridiagonal thomas algorithm
	//until then, Gauss-Seidel
	'Backward Euler + Gauss Seidel' : {
		initOldQ : function() {
			if (this.oldQ == undefined) {
				this.oldQ = [];
				for (var i = 0; i < this.nx; ++i) {
					this.oldQ[i] = [0, 0, 0];
				}
			}
			for (var i = 0; i < this.nx; ++i) {
				for (var j = 0; j < 3; ++j) {
					this.oldQ[i][j] = this.q[i][j];
				}
			}
		},
		initPressure : function() {
			integrationMethods['Backward Euler + Gauss Seidel'].initOldQ.call(this);
		},
		applyMomentumDiffusion : function(dt) {
			for (var iter = 0; iter < gaussSeidelIterations; ++iter) {
				for (var i = 0; i < this.nx; ++i) {
					var u = this.q[i][1] / this.q[i][0];
					var energyTotal = this.q[i][2] / this.q[i][0];
					var energyKinematic = .5 * u * u;
					var energyThermal = energyTotal - energyKinematic;
					this.pressure[i] = (this.gamma - 1) * this.q[i][0] * energyThermal;
				}
				
				for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
					this.q[i][1] = this.oldQ[i][1] - dt * (this.pressure[i+1] - this.pressure[i-1]) / (this.x[i+1] - this.x[i-1]);
				}
			}
		},
		applyWorkDiffusion : function(dt) {
			for (var iter = 0; iter < gaussSeidelIterations; ++iter) {
				for (var i = 0; i < this.nx; ++i) {
					var u = this.q[i][1] / this.q[i][0];
					var energyTotal = this.q[i][2] / this.q[i][0];
					var energyKinematic = .5 * u * u;
					var energyThermal = energyTotal - energyKinematic;
					this.pressure[i] = (this.gamma - 1) * this.q[i][0] * energyThermal;
				}
				
				for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
					var u_inext = this.q[i+1][1] / this.q[i+1][0];
					var u_iprev = this.q[i-1][1] / this.q[i-1][0];
					this.q[i][2] = this.oldQ[i][2] - dt * (this.pressure[i+1] * u_inext - this.pressure[i-1] * u_iprev) / (this.x[i+1] - this.x[i-1]);
				}
			}
		},
		advect : {
			//Linearizing our relationship between current and next timesteps.
			//Treating the flux limiter and the interface velocity as constants when I could consider them in terms of q's. 
			//I get timesteps of .7, when .3 or so is what the max CFL timestep for explicit was giving me,
			// but still see a lot more oscillations in the system.
			Burgers : function(dt) {
				integrationMethods['Backward Euler + Gauss Seidel'].initOldQ.call(this);
			
				for (var iter = 0; iter < gaussSeidelIterations; ++iter) {
					//get velocity at interfaces from state
					for (var ix = this.nghost-1; ix < this.nx+this.nghost-2; ++ix) {
						this.ui[ix] = .5 * (this.q[ix][1] / this.q[ix][0] + this.q[ix-1][1] / this.q[ix-1][0]);
					}
					this.ui[0] = this.ui[this.nx] = 0;

					//compute flux and advect for each state vector
					for (var j = 0; j < 3; ++j) {
						//r_{i-1/2} flux limiter
						for (var i = this.nghost; i < this.nx+this.nghost-3; ++i) {
							var dq = this.q[i][j] - this.q[i-1][j];
							if (Math.abs(dq) > 0) {
								if (this.ui[i] >= 0) {
									this.r[i][j] = (this.q[i-1][j] - this.q[i-2][j]) / dq;
								} else {
									this.r[i][j] = (this.q[i+1][j] - this.q[i][j]) / dq;
								}
							} else {
								this.r[i][j] = 0;
							}
						}
						this.r[0][j] = this.r[1][j] = this.r[this.nx-1][j] = this.r[this.nx][j] = 0;

						//update cells
						
						for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
							var dx = this.xi[i+1] - this.xi[i];
							
							var dflux_qip = 0;
							var dflux_qi = 0;
							var dflux_qin = 0;
							
							//flux limiter
							var phi = fluxMethods[this.fluxMethod](this.r[i][j]);
							if (this.ui[i] >= 0) {
								dflux_qip -= this.ui[i];
							} else {
								dflux_qi -= this.ui[i];
							}
							dflux_qi -= phi * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
							dflux_qip += phi * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
							
							var phi = fluxMethods[this.fluxMethod](this.r[i+1][j]);
							if (this.ui[i+1] >= 0) {
								dflux_qi += this.ui[i+1];
							} else {
								dflux_qin += this.ui[i+1];
							}
							dflux_qin += phi * .5 * Math.abs(this.ui[i+1]) * (1 - Math.abs(this.ui[i+1] * dt / dx));
							dflux_qi -= phi * .5 * Math.abs(this.ui[i+1]) * (1 - Math.abs(this.ui[i+1] * dt / dx));
							
							this.q[i][j] = (this.oldQ[i][j] - dt / dx * (
								dflux_qip * this.q[i-1][j] 
								+ dflux_qin * this.q[i+1][j])) / (1 + dt / dx * dflux_qi);
						}
					}
				}	
			},
			'Riemann / Roe' : function(dt) {
				//TODO
				integrationMethods['Forward Euler'].advect['Riemann / Roe'].call(this, dt);
			}
		}
	}
};

//called with 'this' the HydroState
var advectMethods = {
	Burgers : {
		initStep : function() {},
		calcCFLTimestep : function() {
			var mindum = undefined;
			for (var i = 0; i < this.nx; ++i) {
				var u = this.q[i][1] / this.q[i][0];
				var energyTotal = this.q[i][2] / this.q[i][0];
				var energyKinematic = .5 * u * u;
				var energyThermal = energyTotal - energyKinematic;
				var speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyThermal);
				var dx = this.xi[i+1] - this.xi[i];
				var dum = dx / (speedOfSound + Math.abs(u));
				if (mindum === undefined || dum < mindum) mindum = dum;
			}
			//if (mindum != mindum) throw 'nan';
			return this.cfl * mindum;
		}
	},
	'Riemann / Roe' : {
		initStep : function() {
			//qi[ix] = q_{i-1/2} lies between q_{i-1} = q[i-1] and q_i = q[i]
			//(i.e. qi[ix] is between q[ix-1] and q[ix])
			for (var ix = 1; ix < this.nx; ++ix) {
				//compute Roe averaged interface values
				var densityL = this.q[ix-1][0];
				var velocityL = this.q[ix-1][1] / densityL;
				var energyTotalL = this.q[ix-1][2] / densityL;
				var energyKinematicL = .5 * velocityL * velocityL;
				var energyThermalL = energyTotalL - energyKinematicL;
				var pressureL = (this.gamma - 1) * densityL * energyThermalL;
				var speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
				var hTotalL = energyTotalL + pressureL / densityL;
				var roeWeightL = Math.sqrt(densityL);
				
				var densityR = this.q[ix][0];
				var velocityR = this.q[ix][1] / densityR;
				var energyTotalR = this.q[ix][2] / densityR;
				var energyKinematicR = .5 * velocityR * velocityR;
				var energyThermalR = energyTotalR - energyKinematicR;
				var pressureR = (this.gamma - 1) * densityR * energyThermalR;
				var speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
				var hTotalR = energyTotalR + pressureR / densityR;
				var roeWeightR = Math.sqrt(densityR);
				
				var denom = roeWeightL + roeWeightR;
				var velocity = (roeWeightL * velocityL + roeWeightR * velocityR) / denom;
				var hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;

				//compute eigenvectors and values at the interface based on Roe averages
				buildEigenstate(
					this.interfaceMatrix[ix],
					this.interfaceEigenvalues[ix], 
					this.interfaceEigenvectors[ix], 
					this.interfaceEigenvectorsInverse[ix], 
					velocity, hTotal, this.gamma);
			}
		
		},
		/*
		store eigenvalues and eigenvectors of interfaces
		use the lambdas to calc the DT based on CFL
		*/
		calcCFLTimestep : function() {
	
			var mindum = undefined;
			for (var i = 1; i < this.nx; ++i) {
				var maxLambda = Math.max(0, this.interfaceEigenvalues[i][0], this.interfaceEigenvalues[i][1], this.interfaceEigenvalues[i][2]);
				var minLambda = Math.min(0, this.interfaceEigenvalues[i+1][0], this.interfaceEigenvalues[i+1][1], this.interfaceEigenvalues[i+1][2]);
				var dum = (this.xi[i+1] - this.xi[i]) / (maxLambda - minLambda);
				if (mindum === undefined || dum < mindum) mindum = dum;
			}
			//if (mindum != mindum) throw 'nan';
			return this.cfl * mindum;
		}
	}
};

var HydroState = makeClass({ 
	init : function(args) {
		this.nx = args.size;
		this.gamma = args.gamma;
		this.cfl =.5;
		var x0 = 0;
		var x1 = 100;
		
		//x_i: cell positions
		this.x = new Float32Array(this.nx);
		for (var i = 0; i < this.nx; ++i) {
			this.x[i] = x0 + (x1 - x0) * i / (this.nx-1);
		}
		
		//x_{i-1/2}: interface positions
		this.xi = new Float32Array(this.nx+1);
		for (var i = 1; i < this.nx+1; ++i) {
			this.xi[i] = .5*(this.x[i] + this.x[i-1]);
		}
		this.xi[0] = 2 * this.xi[1] - this.xi[2];
		this.xi[this.nx] = 2 * this.xi[this.nx-1] - this.xi[this.nx-2]; 

		//q_j,i: state vector, stored as q[j][i]
		//q_0,i: density: rho
		//q_1,i: momentum: rho * v
		//q_2,i: work: rho * e
		this.q = [];
		for (var i = 0; i < this.nx; ++i) {
			this.q[i] = [];
		}

		this.resetSod();
		
		//p_i: pressure
		this.pressure = new Float32Array(this.nx);
	
		
		//used for Burgers
		
		
		//r_{i-1/2}	
		this.r = [];
		for (var i = 0; i < this.nx+1; ++i) {
			this.r[i] = [0,0,0];
		}
		
		//f_{i-1/2}: cell flux
		this.flux = [];
		for (var i = 0; i < this.nx+1; ++i) {
			this.flux[i] = [0,0,0];
		}
		
		//only used with Burger's eqn advection code
		//u_{i-1/2}: interface velocity
		this.ui = new Float32Array(this.nx+1);

		//only used with Riemann eqn advection code:
		//calculated before dt
		// used for dt and for fluxTilde
		this.interfaceMatrix = [];
		this.interfaceEigenvalues = [];	//lambda_{i-1/2},j: interface eigenvalues
		this.interfaceEigenvectors = [];		//e_{i-1/2},j,k: interface eigenvectors
		this.interfaceEigenvectorsInverse = [];		//[e_{i-1/2},j,k]^-1: interface eigenvector column matrix inverse 
		for (var ix = 0; ix < this.nx+1; ++ix) {
			this.interfaceMatrix[ix] = [[1,0,0], [0,1,0], [0,0,1]];
			this.interfaceEigenvalues[ix] = [0, 0, 0];
			this.interfaceEigenvectors[ix] = [[1,0,0], [0,1,0], [0,0,1]];
			this.interfaceEigenvectorsInverse[ix] = [[1,0,0], [0,1,0], [0,0,1]];
		}

		//used for Riemann
		this.rTilde = [];
		for (var i = 0; i < this.nx+1; ++i) {
			this.rTilde[i] = [0,0,0];
		}

		//state change over interface in interface Roe eigenspace
		//used for Riemann
		//tilde means in basis of eigenvectors
		this.interfaceDeltaQTilde = [];
		for (var i = 0; i <= this.nx; ++i) {
			this.interfaceDeltaQTilde[i] = [0,0,0];
		}

		//number of ghost cells
		this.nghost = 2;

		//solver configuration
		this.boundaryMethod = 'mirror';
		this.fluxMethod = 'superbee';
		this.advectMethod = 'Riemann / Roe';
		this.integrationMethod = 'Forward Euler';
	},
	resetSod : function() {
		for (var i = 0; i < this.nx; ++i) {
			var x = this.x[i];
			var rho = (x < 30) ? 1 : .1;
			var u = 0;
			var eTotal = 1;
			this.q[i][0] = rho; 
			this.q[i][1] = rho * u; 
			this.q[i][2] = rho * eTotal; 
		}
	},
	resetWave : function() {
		var x0 = this.x[0];
		var x1 = this.x[this.nx-1];
		var xmid = .5 * (x0 + x1);
		var dg = .1 * (x1 - x0);
		for (var i = 0; i < this.nx; ++i) {
			var x = this.x[i];
			var dx = x - xmid;
			var rho = 1 + .3 * Math.exp(-(dx*dx)/(dg*dg));
			var u = 0;
			var eTotal = 1;
			this.q[i][0] = rho;
			this.q[i][1] = rho * u;
			this.q[i][2] = rho * eTotal;
		}
	},
	boundary : function() {
		boundaryMethods[this.boundaryMethod](this.nx, this.q);
	},
	step : function(dt) {
		
		//apply boundary conditions
		this.boundary();
	
		//solve
		integrationMethods[this.integrationMethod].advect[this.advectMethod].call(this, dt);
			
		//boundary again
		this.boundary();

		//compute pressure
		integrationMethods[this.integrationMethod].initPressure.call(this);
	
		//apply momentum diffusion = pressure
		integrationMethods[this.integrationMethod].applyMomentumDiffusion.call(this, dt);

		//apply work diffusion = momentum
		integrationMethods[this.integrationMethod].applyWorkDiffusion.call(this, dt);
	
		//last boundary update
		this.boundary();
	},
	update : function() {
		//do any pre-calcCFLTimestep preparation (Roe computes eigenvalues here)
		advectMethods[this.advectMethod].initStep.call(this, dt);
		
		//get timestep
		var dt;
		if (useCFL) {
			dt = advectMethods[this.advectMethod].calcCFLTimestep.call(this);
		} else {
			dt = fixedDT;
		}

		//do the update
		this.step(dt);
	}
});

var Hydro = makeClass({
	init : function() {
		this.state = new HydroState({
			size : 200,
			gamma : 7/5
		});
	
		//geometry
		this.vertexPositions = new Float32Array(4*this.state.nx);
		this.vertexStates = new Float32Array(6*this.state.nx);
	},
	update : function() {
		//todo adm or something
		//update a copy of the grid and its once-refined
		//...and a once-unrefined ... over mergeable cells only?
		//then test for errors and split when needed
		this.state.update();
	
		//update geometry
		var x = this.state.x;
		var q = this.state.q;
		var nx = this.state.nx;
		for (var i = 0; i < nx; ++i) {
			this.vertexPositions[0+4*i] = x[i];
			this.vertexPositions[1+4*i] = q[i][0]*15;
			this.vertexPositions[2+4*i] = x[i];
			this.vertexPositions[3+4*i] = 0;
			this.vertexStates[0+6*i] = 1;
			this.vertexStates[1+6*i] = 1;
			this.vertexStates[2+6*i] = 1;
			this.vertexStates[3+6*i] = q[i][0];
			this.vertexStates[4+6*i] = Math.abs(q[i][1] / q[i][0]);
			this.vertexStates[5+6*i] = q[i][2] / q[i][0];
		}
	}
});

var hydro = new Hydro();

function update() {
	//iterate
	hydro.update();
	waveVtxBuf.updateData(hydro.vertexPositions);
	waveStateBuf.updateData(hydro.vertexStates);
	//draw
	GL.draw();
	requestAnimFrame(update);
}

function onresize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	//factor out aspectratio from fovY (thus making it fovX)
	var aspectRatio = canvas.width / canvas.height;
	GL.view.fovY = .5 * (xmax - xmin) / aspectRatio;
	GL.view.pos[1] = (ymax + ymin) / 2 / aspectRatio;
	GL.resize();
}

function buildSelect(id, key, map) {
	var select = $('#' + id);
	for (var k in map) {
		var option = $('<option>', {text : k});
		option.appendTo(select);
		if (hydro.state[key] == k) {
			option.attr('selected', 'true');
		}
	}
	select.change(function() {
		hydro.state[key] = select.val();
	});
}

$(document).ready(function(){
	panel = $('#panel');	

	$('#reset-sod').click(function(){ hydro.state.resetSod(); });
	$('#reset-wave').click(function(){ hydro.state.resetWave(); });

	buildSelect('boundary', 'boundaryMethod', boundaryMethods);
	buildSelect('flux-limiter', 'fluxMethod', fluxMethods);
	buildSelect('advect-method', 'advectMethod', advectMethods);
	buildSelect('integration-method', 'integrationMethod', integrationMethods);

	$('#timeStepCFLBased').change(function() {
		if (!$(this).is(':checked')) return;
		useCFL = true;
	});
	$('#timeStepCFL').val(hydro.state.cfl);
	$('#timeStepCFL').change(function() {
		var v = Number($(this).val());
		if (v != v) return;
		hydro.state.cfl = v;
	});
	$('#timeStepFixed').change(function() {
		if (!$(this).is(':checked')) return;
		useCFL = false;
	});
	$('#timeStepValue').val(fixedDT);
	$('#timeStepValue').change(function() {
		var v = Number($(this).val());
		if (v != v) return;	//stupid javascript ... convert anything it doesn't understand to NaNs...
		fixedDT = v;
	});

	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()
	
	try {
		gl = GL.init(canvas);
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}

	GL.view.ortho = true;
	GL.view.zNear = -1;
	GL.view.zFar = 1;
	GL.view.pos[0] = (xmax + xmin) / 2;
	GL.view.pos[1] = (ymax + ymin) / 2;
	
	var plainShader = new GL.ShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex.xy, 0., 1.);
	gl_PointSize = 3.;
}
*/}),
		vertexPrecision : 'best',
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
void main() {
	gl_FragColor = color;
}
*/}),
		fragmentPrecision : 'best',
		uniforms : {
			color : [1,1,1,1]
		}
	});

	//make static grid
	var grid = [];
	for (var i = xmin; i < xmax; i += gridstep) {
		grid.push(i);
		grid.push(ymin);
		grid.push(i);
		grid.push(ymax);
	}
	for (var j = ymin; j < ymax; j += gridstep) {
		grid.push(xmin);
		grid.push(j);
		grid.push(xmax);
		grid.push(j);
	}
	new GL.SceneObject({
		mode : gl.LINES,
		attrs : {
			vertex : new GL.ArrayBuffer({data:grid, dim:2})
		},
		shader : plainShader,
		uniforms : {
			color : [.5,.5,.5,1]
		}
	});

	//make grid
	waveVtxBuf = new GL.ArrayBuffer({
		dim : 2,
		data : hydro.vertexPositions,
		usage : gl.DYNAMIC_DRAW
	});
	waveStateBuf = new GL.ArrayBuffer({
		dim : 3,
		data : hydro.vertexStates,
		usage : gl.DYNAMIC_DRAW
	});
	new GL.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		attrs : {
			vertex : waveVtxBuf,
			state : waveStateBuf
		},
		shader : new GL.ShaderProgram({
			vertexCode : mlstr(function(){/*
attribute vec2 vertex;
attribute vec3 state;
varying vec3 statev;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	statev = state;
	gl_Position = projMat * mvMat * vec4(vertex.xy, 0., 1.);
	gl_PointSize = 3.;
}
*/}),
			vertexPrecision : 'best',
			fragmentCode : mlstr(function(){/*
varying vec3 statev;
void main() {
	gl_FragColor = vec4(statev, 1.);;
}
*/}),
			fragmentPrecision : 'best'
		})
	});

	//start it off
	onresize();
	update();
});
