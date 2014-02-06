/*
1D ADM
2D: Burger at least
2D ADM
2D unstructured
*/

var panel;
var canvas;
var xmin = -1;
var xmax = 1; 
var ymin = -1;
var ymax = 2;
var gridstep = .1;
var useCFL = true;
var fixedDT = .2;
var pause = false;
var gaussSeidelIterations = 20;
var plainShader;
var graphShader;
var gridObj;

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

var copyState = function(srcQ, destQ) {
	if (destQ === undefined) destQ = [];
	for (var i = 0; i < srcQ.length; ++i) {
		if (destQ[i] === undefined) destQ[i] = [];
		for (var j = 0; j < 3; ++j) {
			destQ[i][j] = srcQ[i][j];
		}
	}
	return destQ;
};

var explicitMethods = {
	Euler : function(dt, deriv) {
		var dq_dt = copyState(this.q);
		deriv.call(this, dt, dq_dt);
		addMulState(this.q, dq_dt, dt);
	},
	RK2 : function(dt, deriv) {
		var src = copyState(this.q);
		var k1 = copyState(this.q);
		deriv.call(this, dt, k1);
		addMulState(this.q, k1, .5 * dt);
		var k2 = copyState(this.q);
		deriv.call(this, dt, k2);
		copyState(src, this.q);
		addMulState(this.q, k2, dt);
	},
	RK4 : function(dt, deriv) {
		var src = copyState(this.q);
		var k1 = copyState(this.q);
		deriv.call(this, dt, k1);
		addMulState(this.q, k1, .5 * dt);
		var k2 = copyState(this.q);
		deriv.call(this, dt, k2);
		copyState(src, this.q);
		addMulState(this.q, k2, .5 * dt);
		var k3 = copyState(this.q);
		deriv.call(this, dt, k3);
		copyState(src, this.q);
		addMulState(this.q, k3, dt);
		var k4 = copyState(this.q);
		deriv.call(this, dt, k4);
		copyState(src, this.q);
		addMulState(this.q, k1, dt / 6);
		addMulState(this.q, k2, dt / 3);
		addMulState(this.q, k3, dt / 3);
		addMulState(this.q, k4, dt / 6);
	},
	ICN3 : function(dt, deriv) {
		//first iteration
		var srcQ = copyState(this.q);
		var firstK = copyState(this.q);
		deriv.call(this, dt, firstK);
		addMulState(this.q, firstK, dt);
		var k = copyState(this.q);

		//second and so on
		for (var i = 1; i < 3; ++i) {
			deriv.call(this, dt, k);
			copyState(srcQ, this.q);
			addMulState(this.q, k, .5 * dt);
			addMulState(this.q, firstK, .5 * dt);
		}
	}
};

var addMulState = function(to, from, scalar) {
	for (var i = 0; i < to.length; ++i) {
		for (var j = 0; j < 3; ++j) {
			to[i][j] += scalar * from[i][j];
		}
	}
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
	constant : function(nx,q) {
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
	freeflow : function(nx,q) {
		q[0][0] = q[1][0] = q[2][0];
		q[nx-1][0] = q[nx-2][0] = q[nx-3][0];
		q[0][1] = q[1][1] = q[2][1];
		q[nx-1][1] = q[nx-2][1] = q[nx-3][1];
		q[0][2] = q[1][2] = q[2][2];
		q[nx-1][2] = q[nx-2][2] = q[nx-3][2];
	}
};

var eulerGetPrimitives = function(i) {
	var rho = this.q[i][0];
	var u = this.q[i][1] / rho;
	var e = this.q[i][2] / rho;
	return [rho, u, e];
};

var EulerEquationBurgersSolver = makeClass({
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
	},
	getPrimitives : eulerGetPrimitives
});

/*
time derivative + advection component + source term = 0
d/dt (rho) + d/dx(rho u) = 0
d/dt (rho u) + d/dx(rho u&2) + d/dx (P) = 0
d/dt (rho e_total) + d/dx (rho e_total u) + d/dx (P u) = 0
*/
var EulerEquationBurgersForwardEuler = makeClass({
	super : EulerEquationBurgersSolver,
	step : function(dt) {
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
			for (var i = this.nghost; i < this.nx+1-this.nghost; ++i) {
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
	
		// apply pressure source terms
		
		//boundary again
		this.boundary();

		//compute pressure
		for (var i = 0; i < this.nx; ++i) {
			var u = this.q[i][1] / this.q[i][0];
			var energyTotal = this.q[i][2] / this.q[i][0];
			var energyKinematic = .5 * u * u;
			var energyThermal = energyTotal - energyKinematic;
			this.pressure[i] = (this.gamma - 1) * this.q[i][0] * energyThermal;
		}

		//apply momentum diffusion = pressure
		for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
			this.q[i][1] -= dt * (this.pressure[i+1] - this.pressure[i-1]) / (this.x[i+1] - this.x[i-1]);
		}

		//apply work diffusion = momentum
		for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
			var u_inext = this.q[i+1][1] / this.q[i+1][0];
			var u_iprev = this.q[i-1][1] / this.q[i-1][0];
			this.q[i][2] -= dt * (this.pressure[i+1] * u_inext - this.pressure[i-1] * u_iprev) / (this.x[i+1] - this.x[i-1]);
		}
	}
});

/*
TODO 3x3 block tridiagonal thomas algorithm
until then, Gauss-Seidel

Linearizing our relationship between current and next timesteps.
Treating the flux limiter and the interface velocity as constants when I could consider them in terms of q's. 
I get timesteps of .7, when .3 or so is what the max CFL timestep for explicit was giving me,
 but still see a lot more oscillations in the system.
*/
var EulerEquationBurgersBackwardEulerGaussSeidel = makeClass({
	super : EulerEquationBurgersSolver,
	step : function(dt) {
		for (var i = 0; i < this.nx; ++i) {
			for (var j = 0; j < 3; ++j) {
				this.oldQ[i][j] = this.q[i][j];
			}
		}
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
	
		// only needed for Burgers	
		
		//boundary again
		this.boundary();

		//apply momentum diffusion = pressure
		for (var i = 0; i < this.nx; ++i) {
			for (var j = 0; j < 3; ++j) {
				this.oldQ[i][j] = this.q[i][j];
			}
		}
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

		//apply work diffusion = momentum
		for (var i = 0; i < this.nx; ++i) {
			for (var j = 0; j < 3; ++j) {
				this.oldQ[i][j] = this.q[i][j];
			}
		}
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
	}
});

var EulerEquationBurgersBackwardEulerTridiagonal = makeClass({
	super : EulerEquationBurgersSolver,
	step : function(dt) {
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
				
				this.q[i][j] = (this.q[i][j] - dt / dx * (
					dflux_qip * this.q[i-1][j] 
					+ dflux_qin * this.q[i+1][j])) / (1 + dt / dx * dflux_qi);
			}
		}

		// only needed for Burgers	
		
		//boundary again
		this.boundary();

		var a = [];	//lower band
		var b = [];	//diagonal
		var c = [];	//upper band

		//apply momentum diffusion = pressure
		for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
			var dtOverTwoDx = dt / (this.x[i+1] - this.x[i-1]);
			this.q[i][1] = this.q[i][1] 
				- dtOverTwoDx * (this.gamma - 1) * (this.q[i+1][2] - .5 * this.q[i+1][1] * this.q[i+1][1] / this.q[i+1][0]) 
				+ dtOverTwoDx * (this.gamma - 1) * (this.q[i-1][2] - .5 * this.q[i-1][1] * this.q[i-1][1] / this.q[i-1][0]);
		}

		//apply work diffusion = momentum
		for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
			var dtOverTwoDx = dt / (this.x[i+1] - this.x[i-1]);
			this.q[i][2] = this.q[i][2] 
				- dtOverTwoDx * this.pressure[i+1] * this.q[i+1][1] / this.q[i+1][0]
				+ dtOverTwoDx * this.pressure[i-1] * this.q[i-1][1] / this.q[i-1][0];
		}
	}
});


var GodunovSolver = makeClass({
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
	},
	step : function(dt) {
		var deriv = GodunovSolver.prototype.calcDerivative;
		explicitMethods.RK4.call(this, dt, deriv);
	},
	calcDerivative : function(dt, dq_dt) {
		for (var ix = 1; ix < this.nx; ++ix) {
			for (var j = 0; j < 3; ++j) {
				//the change in state represented in interface eigenbasis
				this.interfaceDeltaQTilde[ix][j] = 
					this.interfaceEigenvectorsInverse[ix][0][j] * (this.q[ix][0] - this.q[ix-1][0])
					+ this.interfaceEigenvectorsInverse[ix][1][j] * (this.q[ix][1] - this.q[ix-1][1])
					+ this.interfaceEigenvectorsInverse[ix][2][j] * (this.q[ix][2] - this.q[ix-1][2]);
			}
		}
		
		for (var j = 0; j < 3; ++j) {
			this.interfaceDeltaQTilde[0][j] = 0;
			this.interfaceDeltaQTilde[this.nx][j] = 0;
		}
		
		for (var ix = this.nghost; ix < this.nx-this.nghost+1; ++ix) {
			for (var j = 0; j < 3; ++j) {
				var interfaceDeltaQTilde = this.interfaceDeltaQTilde[ix][j];
				if (Math.abs(interfaceDeltaQTilde) > 0) {
					if (this.interfaceEigenvalues[ix][j] >= 0) {
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
		var interfaceQL = [];
		var interfaceQR = [];

		var qTildeMid = [];

		//qi[ix] = q_{i-1/2} lies between q_{i-1} = q[i-1] and q_i = q[i]
		//(i.e. qi[ix] is between q[ix-1] and q[ix])
		//Looks good according to "Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro, p.191
		for (var ix = this.nghost-1; ix < this.nx+this.nghost-2; ++ix) {
			//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
			//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
				
			//most sources seem to refer to qL and qR as q_{i-1/2,j,L} = q_{i-1,j} and q{i-1/2,j,R} = q_{i,j}
			//then there is http://bh0.physics.ubc.ca/People/benjamin/projects/nr/report2/node21.html and Hydrodynamics II section 7.5 which say to extrapolate qL and qR based on slope limited stuff 
			
			//trying to go by this: 
			//http://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/A141984.pdf
			//var qL2 = this.q[ix-2];
			var qL = this.q[ix-1];
			var qR = this.q[ix];
			//var qR2 = this.q[ix+1];
		
			//which basii should be used for reprojection?
			//var qTildeL2 = [];
			//var qTildeL = [];
			//var qTildeR = [];
			//var qTildeR2 = [];
			//for (var j = 0; j < 3; ++j) {
				//qTildeL2[j] = 
				//	this.interfaceEigenvectorsInverse[ix-1][0][j] * qL2[0]
				//	+ this.interfaceEigenvectorsInverse[ix-1][1][j] * qL2[1]
				//	+ this.interfaceEigenvectorsInverse[ix-1][2][j] * qL2[2];
				//qTildeL[j] = 
				//	this.interfaceEigenvectorsInverse[ix][0][j] * qL[0]
				//	+ this.interfaceEigenvectorsInverse[ix][1][j] * qL[1]
				//	+ this.interfaceEigenvectorsInverse[ix][2][j] * qL[2];
				//qTildeR[j] = 
				//	this.interfaceEigenvectorsInverse[ix][0][j] * qR[0]
				//	+ this.interfaceEigenvectorsInverse[ix][1][j] * qR[1]
				//	+ this.interfaceEigenvectorsInverse[ix][2][j] * qR[2];
				//qTildeR2[j] = 
				//	this.interfaceEigenvectorsInverse[ix+1][0][j] * qR2[0]
				//	+ this.interfaceEigenvectorsInverse[ix+1][1][j] * qR2[1]
				//	+ this.interfaceEigenvectorsInverse[ix+1][2][j] * qR2[2];
			//}
			
			//q_{i-1/2},j
			//qTildeMid[ix] = [];
			//for (var j = 0; j < 3; ++j) {
			//	qTildeMid[ix][j] = (7 * (qTildeR[j] + qTildeL[j]) - (qTildeL2[j] + qTildeR2[j])) / 12;
			//}
		//}
		//set limits on parabola left and right sides
		//for (var i = 2; i < this.nx-2; ++i) {
			//for (var j = 0; j < 3; ++j) {
				//var qTilde = 
				//var qTildeL = qTildeMid[i][j];
				//var qTildeR = qTildeMid[i+1][j];
			
				//if (a_R,j - a,j) * (a,j - a_L,j) <= 0 then ...
				//  but how do we find a,j if we are using qTilde i.e. eigenbasis-projected values?
				//  should we use the eigenbasis at the cell center?  or the eigenbasis at the interfaces somehow?
				//...or should we not be doing this in eigenbasis space?
				//if (qTildeR - 
			//}
		//}
		//for (var ix = 1; ix < this.nx; ++ix) {
			//... and then apply the rest of the interface jacobian 
			//for (var j = 0; j < 3; ++j) {
			//	interfaceQL[j] = 
			//		this.interfaceEigenvectors[ix][0][j] * this.interfaceEigenvalues[ix][0] * qTildeL[0]
			//		+ this.interfaceEigenvectors[ix][1][j] * this.interfaceEigenvalues[ix][1] * qTildeL[1]
			//		+ this.interfaceEigenvectors[ix][2][j] * this.interfaceEigenvalues[ix][2] * qTildeL[2];
			//	interfaceQR[j] = 
			//		this.interfaceEigenvectors[ix][0][j] * this.interfaceEigenvalues[ix][0] * qTildeR[0]
			//		+ this.interfaceEigenvectors[ix][1][j] * this.interfaceEigenvalues[ix][1] * qTildeR[1]
			//		+ this.interfaceEigenvectors[ix][2][j] * this.interfaceEigenvalues[ix][2] * qTildeR[2];
			//	fluxAvg[j] = .5 * (interfaceQL[j] + interfaceQR[j]);
			//}
		//}
		//for (var ix = 1; ix < this.nx; ++ix) {

/**/
			//if I wasn't doing all the above slope-limited stuff, this would suffice:
			//however if I'm not using this then I don't need to store the interface jacobian matrix
			for (var j = 0; j < 3; ++j) {
				fluxAvg[j] = .5 * ( 
					this.interfaceMatrix[ix][0][j] * (qL[0] + qR[0])
					+ this.interfaceMatrix[ix][1][j] * (qL[1] + qR[1])
					+ this.interfaceMatrix[ix][2][j] * (qL[2] + qR[2]));
			}
/**/

			//calculate flux
			for (var j = 0; j < 3; ++j) {
				var theta = 0;
				if (this.interfaceEigenvalues[ix][j] >= 0) {
					theta = 1;
				} else {
					theta = -1;
				}
				
				var phiTilde = fluxMethods[this.fluxMethod](this.rTilde[ix][j]);
				var dx = this.xi[ix] - this.xi[ix-1];
				var epsilon = this.interfaceEigenvalues[ix][j] * dt / dx;
				
				//flux[ix][k] = fluxTilde[ix][j] * interfaceEigenvectors[ix][k][j]
				//flux in eigenvector basis is the q vector transformed by the inverse then scaled by the eigenvalue
				//should the eigenvalue be incorperated here, after flux limiter is taken into account, or beforehand?
				//1D says after, but notes say before ...
				var deltaFluxTilde = this.interfaceEigenvalues[ix][j] * this.interfaceDeltaQTilde[ix][j];
				
				fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phiTilde * (epsilon - theta));
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
		if (dq_dt === undefined) dq_dt = [];
		for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
			if (dq_dt[i] === undefined) dq_dt[i] = [];
			for (var j = 0; j < 3; ++j) {
				//this.q[i][j] -= dt * (this.flux[i+1][j] - this.flux[i][j]) / (this.xi[i+1] - this.xi[i]);
				dq_dt[i][j] = -(this.flux[i+1][j] - this.flux[i][j]) / (this.xi[i+1] - this.xi[i]);
			}
		}
		for (var i = 0; i < this.nghost; ++i) {
			dq_dt[i] = [0,0,0];
			dq_dt[this.nx-i-1] = [0,0,0];
		}
	}
});

var EulerEquationGodunovSolver = makeClass({
	super : GodunovSolver,
	buildEigenstate : function(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocity, hTotal, gamma) {
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
	},
	getPrimitives : eulerGetPrimitives
});

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
var EulerEquationGodunovForwardEuler = makeClass({
	super : EulerEquationGodunovSolver,
	initStep : function() {
		//qi[ix] = q_{i-1/2} lies between q_{i-1} = q[i-1] and q_i = q[i]
		//(i.e. qi[ix] is between q[ix-1] and q[ix])
		for (var ix = 1; ix < this.nx; ++ix) {
			//compute averaged interface values
			var densityL = this.q[ix-1][0];
			var velocityL = this.q[ix-1][1] / densityL;
			var energyTotalL = this.q[ix-1][2] / densityL;
			var energyKinematicL = .5 * velocityL * velocityL;
			var energyThermalL = energyTotalL - energyKinematicL;
			var pressureL = (this.gamma - 1) * densityL * energyThermalL;
			var speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
			var hTotalL = energyTotalL + pressureL / densityL;
			
			var densityR = this.q[ix][0];
			var velocityR = this.q[ix][1] / densityR;
			var energyTotalR = this.q[ix][2] / densityR;
			var energyKinematicR = .5 * velocityR * velocityR;
			var energyThermalR = energyTotalR - energyKinematicR;
			var pressureR = (this.gamma - 1) * densityR * energyThermalR;
			var speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
			var hTotalR = energyTotalR + pressureR / densityR;
			
			var velocity = .5 * (velocityL + velocityR);
			var hTotal = .5 * (hTotalL + hTotalR);

			//compute eigenvectors and values at the interface based on averages
			EulerEquationGodunovForwardEuler.prototype.buildEigenstate.call(this,
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix], 
				velocity, hTotal, this.gamma);
		}	
	}
});

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
var EulerEquationRoeForwardEuler = makeClass({
	super : EulerEquationGodunovSolver,
	initStep : function() {
		//same idea as Godunov but with Roe weighting: sqrt(rho)
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
			EulerEquationRoeForwardEuler.prototype.buildEigenstate.call(this,
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix], 
				velocity, hTotal, this.gamma);
		}
	}
});

var eulerEquationSimulation = {
	methods : {
		'Burgers / Forward Euler' : EulerEquationBurgersForwardEuler.prototype,
		'Burgers / Backward Euler via Gauss Seidel' : EulerEquationBurgersBackwardEulerGaussSeidel.prototype,
		//'Burgers / Backward Euler via Block Tridiagonal' : EulerEquationBurgersBackwardEulerTridiagonal.prototype,
		'Godunov / Forward Euler' : EulerEquationGodunovForwardEuler.prototype,
		'Roe / Forward Euler' : EulerEquationRoeForwardEuler.prototype
	},
	initialConditions : {
		Sod : function() {
			this.resetCoordinates(-1, 1);
			for (var i = 0; i < this.nx; ++i) {
				var x = this.x[i];
				var rho = (x < (xmin * .7 + xmax * .3)) ? 1 : .1;
				var u = 0;
				var eTotal = 1;
				this.q[i][0] = rho; 
				this.q[i][1] = rho * u; 
				this.q[i][2] = rho * eTotal; 
			}
		},
		Sedov : function() {
			this.resetCoordinates(-1, 1);
			var pressure = 1e-5;
			var rho = 1;
			for (var i = 0; i < this.nx; ++i) {
				this.q[i][0] = rho;
				this.q[i][1] = 0;
				this.q[i][2] = pressure / (this.gamma - 1);
			}
			this.q[Math.floor(this.nx/2)][2] = 1e+5;
		},
		Advect : function() {
			this.resetCoordinates(-1, 1);
			var xmid = .5 * (xmax + xmin);
			for (var i = 0; i < this.nx; ++i) {
				var x = this.x[i];
				var xGreaterThanMid = x > xmid;
				var rho = xGreaterThanMid ? .5 : 1;
				var u = 1;
				var energyKinetic = .5 * u * u;
				var pressure = 1;
				var energyTotal = pressure / (rho * (this.gamma - 1)) + energyKinetic;
				this.q[i][0] = rho;
				this.q[i][1] = rho * u;
				this.q[i][2] = rho * energyTotal;
			}
		},
		Wave : function() {
			this.resetCoordinates(-1, 1);
			var xmid = .5 * (xmin + xmax);
			var dg = .1 * (xmax - xmin);
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
		}
	}
};

/*
From "Numerical Relativity", Alcubierre. 2008. ch.10 "Examples of Numerical Relativity" section 2: "Toy 1+1 Relativity"

State variables:
q[0] = D_alpha = d/dx ln alpha
q[1] = D_g = d/dx ln g
q[2] = K^tilde = sqrt(g) K

... for primtiive variables
alpha = lapse
g = g_xx = metric
K = extrinsic curvature = K^x_x

...and some extra variables:
f comes from Bona-Masso slicing family: d/dt alpha = -alpha^2 f(alpha) K
(shift value is zero)
*/

var adm_BonaMasso_f = 1; 		//aiming to recreate the f=.5, f=1, f=1.5 graphs ...

function admEquationsBuildEigenstate(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, alpha, g) {
	var f = adm_BonaMasso_f;
	
	var oneOverSqrtG = 1 / Math.sqrt(g);
	var sqrtF = Math.sqrt(f);
	
	matrix[0][0] = 0;
	matrix[0][1] = 0;
	matrix[0][2] = alpha * oneOverSqrtG;
	matrix[1][0] = 0;
	matrix[1][1] = 0;
	matrix[1][2] = 0;
	matrix[2][0] = alpha * f * oneOverSqrtG;
	matrix[2][1] = 2 * alpha * oneOverSqrtG;
	matrix[2][2] = 0;

	eigenvalues[0] = -alpha * sqrtF * oneOverSqrtG;
	eigenvalues[1] = 0;
	eigenvalues[2] = -eigenvalues[0];

	//column 0: min
	eigenvectors[0][0] = f;
	eigenvectors[0][1] = 2;
	eigenvectors[0][2] = -sqrtF;

	//column 1: mid
	eigenvectors[1][0] = 0;
	eigenvectors[1][1] = 1;
	eigenvectors[1][2] = 0;

	//column 2: max
	eigenvectors[2][0] = f;
	eigenvectors[2][1] = 2;
	eigenvectors[2][2] = sqrtF;
	
	//calculate eigenvector inverses ... 
	mat33invert(eigenvectorsInverse, eigenvectors);
}

var ADMEquationGodunovForwardEuler = makeClass({
	super : GodunovSolver,
	initStep : function() {
		var dx = (xmax - xmin) / this.nx;
		//same idea as Godunov but with Roe weighting: sqrt(rho)
		for (var ix = 2; ix < this.nx-1; ++ix) {
			var ixL = ix-1;
			var ixR = ix;

			//q_ix,0 = d/dx ln alpha
			var ln_alpha_L = (this.q[ixL+1][0] - this.q[ixL-1][0]) / (2 * dx);
			var alpha_L = Math.exp(ln_alpha_L);

			//q_ix,1 = d/dx ln g
			var ln_g_L = (this.q[ixL+1][0] - this.q[ixL-1][0]) / (2 * dx);
			var g_L = Math.exp(ln_g_L);

			var ln_alpha_R = (this.q[ixR+1][0] - this.q[ixR-1][0]) / (2 * dx);
			var alpha_R = Math.exp(ln_alpha_R);
			
			var ln_g_R = (this.q[ixR+1][0] - this.q[ixR-1][0]) / (2 * dx);
			var g_R = Math.exp(ln_g_R);

			var alpha = .5 * (alpha_L + alpha_R);
			var g = .5 * (g_L + g_R);

			//compute eigenvectors and values at the interface based on Roe averages
			admEquationsBuildEigenstate(
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix], 
				alpha, g);
		}
		//how about those boundary eigenstates?
	},
	step : function(dt) {
		var deriv = ADMEquationGodunovForwardEuler.prototype.calcDerivative;

		//ICN3:
		//first iteration
		var srcQ = copyState(this.q);
		var firstK = deriv.call(this, dt);
		addMulState(this.q, firstK, dt);

		//second and so on
		for (var i = 1; i < 3; ++i) {
			var k = deriv.call(this, dt);
			copyState(srcQ, this.q);
			addMulState(this.q, k, .5 * dt);
			addMulState(this.q, firstK, .5 * dt);
		}
	},
	getPrimitives : function(i) {
		return this.q[i];
	}
});

var admEquationSimulation = {
	methods : {
		'Godunov / Forward Euler' : ADMEquationGodunovForwardEuler.prototype
	},
	initialConditions : {
		GaugeShock : function() {
			this.resetCoordinates(-30, 30);
			var xmid = (xmax + xmin) * .5;
			for (var i = 0; i < this.nx; ++i) {
				var x = (this.x[i] - xmid) / ((xmax - xmid) / 3);
				var h = Math.exp(-x*x); 
				var dh_dx = -2 * x * h;
				var d2h_dx2 = 2 * h * (2 * x * x - 1);
				var g = 1 - dh_dx * dh_dx;
				var D_g = -2 * dh_dx * d2h_dx2 / g;
				var KTilde = -d2h_dx2 / g;
				var f = adm_BonaMasso_f;
				var D_alpha = Math.sqrt(f) * KTilde;
				this.q[i][0] = D_alpha;
				this.q[i][1] = D_g;
				this.q[i][2] = KTilde;
			}
		}
	}
};

//hmm, maybe I should combine all of these into one list, and have them individuall state who they belong to
var simulations = {
	Euler : eulerEquationSimulation,
	ADM : admEquationSimulation
};

var HydroState = makeClass({ 
	//call this to reflect this.x and this.xi after changing xmin, xmax, or this.nx
	resetCoordinates : function(xmin_, xmax_) {
		xmin = xmin_;
		xmax = xmax_;
		onresize();

		//x_i: cell positions
		this.x = new Float32Array(this.nx);
		for (var i = 0; i < this.nx; ++i) {
			this.x[i] = xmin + (xmax - xmin) * i / (this.nx-1);
		}
		
		//x_{i-1/2}: interface positions
		this.xi = new Float32Array(this.nx+1);
		for (var i = 1; i < this.nx+1; ++i) {
			this.xi[i] = .5*(this.x[i] + this.x[i-1]);
		}
		this.xi[0] = 2 * this.xi[1] - this.xi[2];
		this.xi[this.nx] = 2 * this.xi[this.nx-1] - this.xi[this.nx-2]; 
	},
	init : function(args) {
		this.nx = args.size;
		this.cfl =.5;
		
		//used for Euler equation of state
		this.gamma = args.gamma;
	
		this.resetCoordinates();

		//q_j,i: state vector, stored as q[j][i]
		//q_0,i: density: rho
		//q_1,i: momentum: rho * v
		//q_2,i: work: rho * e
		this.q = [];
		for (var i = 0; i < this.nx; ++i) {
			this.q[i] = [];
		}

		eulerEquationSimulation.initialConditions.Sod.call(this);
		
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

	
		//used for Backward Euler + Gauss Seidel
		this.oldQ = [];
		for (var i = 0; i < this.nx; ++i) {
			this.oldQ[i] = [0, 0, 0];
		}
		
		//number of ghost cells
		this.nghost = 2;

		//solver configuration
		this.boundaryMethod = 'mirror';
		this.fluxMethod = 'superbee';
		this.simulation = 'Euler';
		this.algorithm = 'Roe / Forward Euler';
		this.explicitMethod = 'RK4';
	},
	boundary : function() {
		boundaryMethods[this.boundaryMethod](this.nx, this.q);
	},
	step : function(dt) {
		//apply boundary conditions
		this.boundary();
	
		//solve
		simulations[this.simulation].methods[this.algorithm].step.call(this, dt);
	},
	getPrimitives : function() {
		return simulations[this.simulation].methods[this.algorithm].getPrimitives.apply(this, arguments);
	},
	update : function() {
		//do any pre-calcCFLTimestep preparation (Roe computes eigenvalues here)
		simulations[this.simulation].methods[this.algorithm].initStep.call(this, dt);
		
		//get timestep
		var dt;
		if (useCFL) {
			dt = simulations[this.simulation].methods[this.algorithm].calcCFLTimestep.call(this);
		} else {
			dt = fixedDT;
		}
window.lastDT = dt;

		//do the update
		this.step(dt);
	}
});

var Hydro = makeClass({
	init : function() {
		var size = Number($.url().param('size'));
		if (size === undefined || size !== size) size = 200;
		var gamma = Number($.url().param('gamma'));
		if (gamma === undefined || gamma !== gamma) gamma = 7/5;

		this.state = new HydroState({
			size : size,
			gamma : gamma
		});
		
		//geometry buffers
		this.vertexXBuffer = new GL.ArrayBuffer({
			dim : 1,
			count : size,
			usage : gl.DYNAMIC_DRAW,
			keep : true
		});

		var colors = [[1,0,0,1], [0,1,0,1], [0,.5,1,1], [1,1,1,1]];

		this.primitiveBuffers = [];
		this.stateGraphObjs = [];
		for (var i = 0; i < 3; ++i) {
			var stateBuffer = new GL.ArrayBuffer({
				dim : 1,
				count : size,
				usage : gl.DYNAMIC_DRAW,
				keep : true
			});
			this.primitiveBuffers.push(stateBuffer);
	
			//make graphs
			var stateGraphObj = new GL.SceneObject({
				mode : gl.LINE_STRIP,
				attrs : {
					vertex : this.vertexXBuffer,
					state : stateBuffer
				},
				shader : graphShader,
				uniforms : {
					color : colors[i % colors.length]
				},
				blend : [gl.ONE, gl.ONE]
			});
			this.stateGraphObjs.push(stateGraphObj);
		}
	},
	update : function() {
		//update a copy of the grid and its once-refined
		//...and a once-unrefined ... over mergeable cells only?
		//then test for errors and split when needed
		this.state.update();
	
		//update geometry
		for (var i = 0; i < this.state.nx; ++i) {
			this.vertexXBuffer.data[i] = this.state.x[i];
			//rescale to -1,1
			this.vertexXBuffer.data[i] -= xmin;
			this.vertexXBuffer.data[i] *= 2 / (xmax - xmin);
			this.vertexXBuffer.data[i]--;
		}
		this.vertexXBuffer.updateData();
	
		for (var i = 0; i < this.stateGraphObjs.length; ++i) {
			this.stateGraphObjs[i].uniforms.scale = 1/(xmax - xmin);
		}

		//TODO this is only for Euler -- make it generic
		for (var i = 0; i < this.state.nx; ++i) {
			var prims = this.state.getPrimitives(i);
			this.primitiveBuffers[0].data[i] = prims[0];
			this.primitiveBuffers[1].data[i] = prims[1];
			this.primitiveBuffers[2].data[i] = prims[2];
		}
		for (var j = 0; j < 3; ++j) {
			this.primitiveBuffers[j].updateData();
		}
	}
});

var hydro;

function update() {
	//iterate
	if (!pause) hydro.update();
	//draw
	GL.draw();
	requestAnimFrame(update);
}

function onresize() {
	GL.canvas.width = window.innerWidth;
	GL.canvas.height = window.innerHeight;
	$('#content').height(window.innerHeight - 50);
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
	var canvas = $('<canvas>', {
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
	GL.view.pos[0] = 0;//(xmax + xmin) / 2;
	GL.view.pos[1] = 0;//(ymax + ymin) / 2;
	GL.view.fovY = ymax - ymin;

	plainShader = new GL.ShaderProgram({
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

	graphShader = new GL.ShaderProgram({
		vertexCode : mlstr(function(){/*
attribute float vertex;
attribute float state;
uniform float scale;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex, scale * state, 0., 1.);
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
			color : [1,1,1,1],
			scale : 1
		}
	});

	hydro = new Hydro();
	
	panel = $('#panel');	

	var panelContent = $('#content');
	$('#menu').click(function() {
		if (panelContent.css('display') == 'none') {
			panelContent.show();
		} else {
			panelContent.hide();
		}
	});

	$('#pause').click(function() {
		pause = !pause;
		$('#pause').attr('src', pause ? 'play.png' : 'pause.png');
	});

	for (simulationName in simulations) {
		var parent = $('#'+simulationName+'_controls');
		var simulation = simulations[simulationName];
		for (initialConditionName in simulation.initialConditions) {
			(function(){
				var method = simulation.initialConditions[initialConditionName];
				$('<button>', {
					text : 'Reset '+initialConditionName,
					click : function() {
						method.call(hydro.state);
					}
				}).appendTo(parent);
				$('<br>').appendTo(parent);
			})();
		}
	}

	(function(){
		var select = $('#gridsize');
		$.each([20, 50, 100, 200, 500, 1000], function(i,gridsize){
			var option = $('<option>', {text : gridsize});
			if (hydro.state.nx == gridsize) option.attr('selected', 'true');
			option.appendTo(select);
		});
		select.change(function(){
			var params = $.url().param();
			params.size = select.val();
			var url = location.href.match('[^?]*');
			var sep = '?';
			for (k in params) {
				if (k == '' && params[k] == '') continue;
				url += sep;
				url += k + '=' + params[k];
				sep = '&';
			}
			location.href = url;
		});
	})();
	
	buildSelect('boundaryMethod', 'boundaryMethod', boundaryMethods);
	buildSelect('fluxMethod', 'fluxMethod', fluxMethods);
	
	buildSelect('Euler_algorithm', 'algorithm', eulerEquationSimulation.methods);	//this will change as 'simulations' changes
	buildSelect('ADM_algorithm', 'algorithm', admEquationSimulation.methods);	//this will change as 'simulations' changes
	
	buildSelect('explicitMethod', 'explicitMethod', explicitMethods);

	(function(){
		var id = 'simulation';
		var key = 'simulation';
		var map = simulations;
		var select = $('#' + id);
		for (var k in map) {
			var option = $('<option>', {text : k});
			option.appendTo(select);
			if (hydro.state[key] == k) {
				option.attr('selected', 'true');
			}
		}
		select.change(function() {
			var val = select.val();
			hydro.state[key] = val; 
			
			$('#ADM_controls, #Euler_controls').hide();
			$('#'+val+'_controls').show();
		
			//...and now select the visible select's selection
			hydro.state.algorithm = $('#'+val+'_algorithm').val();
		});
		
		select.val(hydro.state.simulation);
		select.trigger('change');
	})();

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


	var axisColor = [.75, .75, .75, 1];
	(new GL.SceneObject({
		mode : gl.LINES,
		attrs : {
			vertex : new GL.ArrayBuffer({dim : 2, data : [xmin, 0, xmax, 0]})
		},
		shader : plainShader,
		uniforms : {
			color : axisColor
		}
	})).prependTo(GL.root);

	(new GL.SceneObject({
		mode : gl.LINES,
		attrs : {
			vertex : new GL.ArrayBuffer({dim : 2, data : [0, ymin, 0, ymax]})
		},
		shader : plainShader,
		uniforms : {
			color : axisColor
		}
	})).prependTo(GL.root);


	//make static grid
	//...might not be so static...
	var grid = [];
	for (var i = Math.floor(xmin / gridstep) * gridstep; 
		i <= Math.ceil(xmax / gridstep) * gridstep; 
		i += gridstep) 
	{
		grid.push(i);
		grid.push(ymin);
		grid.push(i);
		grid.push(ymax);
	}
	for (var j = Math.floor(ymin / gridstep) * gridstep; 
		j < Math.ceil(ymax / gridstep) * gridstep; 
		j += gridstep)
	{
		grid.push(xmin);
		grid.push(j);
		grid.push(xmax);
		grid.push(j);
	}
	gridObj = new GL.SceneObject({
		mode : gl.LINES,
		attrs : {
			vertex : new GL.ArrayBuffer({data:grid, dim:2})
		},
		shader : plainShader,
		uniforms : {
			color : [.25,.25,.25,1]
		}
	});
	gridObj.prependTo(GL.root);
	


	var zoomFactor = .0003;	// upon mousewheel
	var dragging = false;
	mouse = new Mouse3D({
		pressObj : canvas,
		mousedown : function() {
			dragging = false;
		},
		move : function(dx,dy) {
			dragging = true;
			var aspectRatio = canvas.width / canvas.height;
			GL.view.pos[0] -= dx / canvas.width * 2 * (aspectRatio * GL.view.fovY);
			GL.view.pos[1] += dy / canvas.height * 2 * GL.view.fovY;
			GL.updateProjection();
		},
		zoom : function(zoomChange) {
			dragging = true;
			var scale = Math.exp(-zoomFactor * zoomChange);
			GL.view.fovY *= scale 
			GL.updateProjection();
		}
	});
	
	//start it off
	onresize();
	$(window).resize(onresize);
	update();
});
