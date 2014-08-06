/*
sources:
http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/
http://www.cfdbooks.com/cfdcodes.html
"Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro
http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf also does not

lots of accuracy issues with the GPU version ... or bugs I'm not finding
in case of accuracy issues, check out view-source:http://hvidtfeldts.net/WebGL-DP/webgl.html for vec2 single -> double encoding 
*/

var gl;
var glutil;
var panel;
var canvas;

var xmin = -.5;
var xmax = .5; 
var ymin = -.5;
var ymax = .5;

var useNoise = true;
var noiseAmplitude = .01;

var useCFL = true;
var fixedDT = .00025;

var mouse;

var externalForceX = 0;
var externalForceY = 0;

var fbo;
var FloatTexture2D;

var quadObj, lineObj;
		
var allShaders = [];

var drawToScreenShader = {};

var solidShader;
var copyShader;
var minReduceShader;

var resetSodShader;
var resetSodCylinderSolidShader;
var resetWaveShader;
var resetKelvinHemholtzShader;

var burgersComputeCFLShader;
var burgersComputeInterfaceVelocityShader;
var burgersComputeFluxSlopeShader = [];	//[side]
var burgersComputeFluxShader = {};		//[fluxMethod][side]
var burgersUpdateStateShader;
var burgersComputePressureShader;
var burgersApplyPressureToMomentumShader;
var burgersApplyPressureToWorkShader;

var roeComputeCFLShader;
var roeComputeRoeValueShader = [];	//[side] = [velocity.x, velocity.y, hTotal, speedOfSound]
var roeComputeEigenvalueShader = [];		//[side]
var roeComputeEigenvectorColumnShader = [];	//[side][column state]
var roeComputeEigenvectorInverseColumnShader = [];	//[side][column state]
var roeComputeDeltaQTildeShader = [];	//[side]
var roeComputeFluxSlopeShader = [];	//[side]
var roeComputeFluxShader = {};	//[fluxMethod][side]
var roeUpdateStateShader;

var encodeTempTex;
var encodeShader = [];	//[channel]

GLUtil.prototype.oninit.push(function() {
	var glutil = this;
	glutil.KernelShader = makeClass({
		super : glutil.ShaderProgram,
		init : function(args) {
			
			var varyingCodePrefix = 'varying vec2 pos;\n';

			var fragmentCodePrefix = '';
			var uniforms = {};
			if (args.uniforms !== undefined) {
				$.each(args.uniforms, function(uniformName, uniformType) {
					if ($.isArray(uniformType)) {
						//save initial value
						uniforms[uniformName] = uniformType[1];
						uniformType = uniformType[0];
					}
					fragmentCodePrefix += 'uniform '+uniformType+' '+uniformName+';\n';
				});
			}
			if (args.texs !== undefined) {
				for (var i = 0; i < args.texs.length; ++i) {
					var v = args.texs[i];
					var name, vartype;
					if (typeof(v) == 'string') {
						name = v;
						vartype = 'sampler2D';
					} else {
						name = v[0];
						vartype = v[1];
					}
					fragmentCodePrefix += 'uniform '+vartype+' '+name+';\n';
					uniforms[name] = i;
				}
			}


			if (!glutil.KernelShader.prototype.kernelVertexShader) {
				glutil.KernelShader.prototype.kernelVertexShader = new glutil.VertexShader({
					code : 
						glutil.vertexPrecision + 
						varyingCodePrefix +
						mlstr(function(){/*
attribute vec2 vertex;
attribute vec2 texCoord;
void main() {
	pos = texCoord; 
	gl_Position = vec4(vertex, 0., 1.);
}
*/})
				});	
			}

			args.vertexShader = glutil.KernelShader.prototype.kernelVertexShader;
			args.fragmentCode = glutil.fragmentPrecision + varyingCodePrefix + fragmentCodePrefix + args.code;
			delete args.code;
			args.uniforms = uniforms;	
			glutil.KernelShader.super.call(this, args);
		}
	});
});

//coordinate names
var coordNames = ['x', 'y'];

var drawToScreenMethods = {
	Density : 'return texture2D(qTex, pos).x;',
	Velocity : 'vec4 q = texture2D(qTex, pos); return length(q.yz) / q.x;',
	Energy : 'vec4 q = texture2D(qTex, pos); return q.w / q.x;',
	//P = (gamma - 1) rho (eTotal - eKinetic)
	Pressure : 'return texture2D(pressureTex, pos).x;',
	//Curl : 'vec4 q = texture2D(qTex, pos); return (dFdy(q.y) * (rangeMax.x - rangeMin.x) * dpos.x - dFdx(q.z) * (rangeMax.y - rangeMin.y) * dpos.y) / q.w;'
	Curl : 'vec4 q = texture2D(qTex, pos); return (dFdy(q.y)  - dFdx(q.z)) / q.w;'
};

var fluxMethods = {
	'donor cell' : 'return vec4(0., 0., 0., 0.);',
	'Lax-Wendroff' : 'return vec4(1., 1., 1., 1.);',
	'Beam-Warming' : 'return r;',
	'Fromm' : 'return .5 * (1. + r);',

	//Wikipedia
	//CHARM : function(r) { return Math.max(0, r*(3*r+1)/((r+1)*(r+1)) ); },
	//HCUS : function(r) { return Math.max(0, 1.5 * (r + Math.abs(r)) / (r + 2) ); },
	//HQUICK : function(r) { return Math.max(0, 2 * (r + Math.abs(r)) / (r + 3) ); },
	//Koren : function(r) { return Math.max(0, Math.min(2*r, (1 + 2*r)/3 ,2) ); },
	//minmod : function(r) { return Math.max(0, Math.min(r,1) ); },
	//Oshker : function(r) { return Math.max(0, Math.min(r,1.5) ); },	//replace 1.5 with 1 <= beta <= 2	
	//ospre : function(r) { return .5 * (r*r + r) / (r*r + r + 1); },
	//smart : function(r) { return Math.max(0, Math.min(2 * r, .25 + .75 * r, 4)); },
	//Sweby : function(r) { return Math.max(0, Math.min(1.5 * r, 1), Math.min(r, 1.5)); },	//replace 1.5 with 1 <= beta <= 2
	//UMIST : function(r) { return Math.max(0, Math.min(2*r, .75 + .25*r, .25 + .75*r, 2)); },	
	//'van Albada 1' : function(r) { return (r * r + r) / (r * r + 1); },
	//'van Albada 2' : function(r) { return 2 * r / (r * r + 1); },
	
	'van Leer' : 'return (r + abs(r)) / (1. + abs(r));', 
	'monotonized central' : 'return max(vec4(0., 0., 0., 0.), min(vec4(2.), min(.5 * (1. + r), 2. * r)));',
	superbee : 'return max(vec4(0., 0., 0., 0.), max(min(vec4(1.), 2. * r), min(vec4(2.), r)));'
	
	//'Barth-Jespersen' : function(r) {.5 * (r + 1) * Math.min(1, 4*r/(r+1), 4/(r+1)); }
};

//'this' is HydroState
function drawBoundaries(color) {
	var nx = this.nx;
	//left
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			.5/nx, 0,
			.5/nx, 1
		]
	});
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			1.5/nx, 0,
			1.5/nx, 1
		]
	});	
	//right
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			(nx-1.5)/nx, 0,
			(nx-1.5)/nx, 1
		]
	});
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			(nx-.5)/nx, 0,
			(nx-.5)/nx, 1
		]
	});
	//bottom
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			0, .5/nx,
			1, .5/nx
		]
	});
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			0, 1.5/nx,
			1, 1.5/nx
		]
	});	
	//top
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			0, (nx-1.5)/nx,
			1, (nx-1.5)/nx
		]
	});
	this.drawLine({
		shader : solidShader,
		uniforms : {
			color : color
		},
		vtxs : [
			0, (nx-.5)/nx,
			1, (nx-.5)/nx
		]
	});
}

var boundaryMethods = {
	//currently an error with this -- either from pixels being left, or a scratch tex bug in general,
	// but the lower-right corner diverges and explodes shortly after switching boundaries on and then off again
	//...and the problem goes away when the drawBoundaries goes away
	periodic : function() {
		//clear borders
		fbo.bind();
		gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.scratchTex.obj, 0);
		fbo.check();
		quadObj.draw({
			shader : copyShader,
			texs : [this.solidTex]
		});
		drawBoundaries.call(this, [0,0,0,0]);
		fbo.unbind();
		this.swapTexs('scratchTex', 'solidTex');
	},
	mirror : function(nx,q) {
		//set borders
		fbo.bind();
		gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.scratchTex.obj, 0);
		fbo.check();
		quadObj.draw({
			shader : copyShader,
			texs : [this.solidTex]
		});
//drawBoundaries.call(this, [0,0,0,0]);
		drawBoundaries.call(this, [1,1,1,1]);
		fbo.unbind();
		this.swapTexs('scratchTex', 'solidTex');
	},
	/*
	constant : function() {
		//TODO set all lookup texture wraps to ... zero?
	},
	*/
	freeflow : function(nx,q) {
		//set all lookup texture wraps to ... CLAMP
		$.each(this.allFloatTexs, function(i, tex) {
			tex.bind();
			tex.setWrap({s : gl.CLAMP_TO_EDGE, t : gl.CLAMP_TO_EDGE});
		});
	}
};

//called with 'this' the HydroState
var advectMethods = {
	Burgers : {
		initStep : function() {
		},
		calcCFLTimestep : function() {
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.cflReduceTex.obj, 0);
			fbo.check();
			quadObj.draw({
				shader : burgersComputeCFLShader,
				uniforms : {
					dpos : [1/this.nx, 1/this.nx],
					rangeMin : [xmin, ymin],
					rangeMax : [xmax, ymax],
					externalForce : [externalForceX, externalForceY],
					gamma : this.gamma
				},
				texs : [this.qTex, this.solidTex]
			});
			fbo.unbind();	
			return this.reduceToDetermineCFL();
		},	
		advect : function(dt) {			
			var dx = (xmax - xmin) / this.nx;
			var dy = (ymax - ymin) / this.nx;
			var dxi = [dx, dy];
			
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.uiTex.obj, 0);
			fbo.check();
			//get velocity at interfaces from state
			quadObj.draw({
				shader : burgersComputeInterfaceVelocityShader,
				uniforms : {
					dpos : [1/this.nx, 1/this.nx]
				},
				texs : [this.qTex, this.solidTex]
			});
			fbo.unbind();
			for (var side = 0; side < 2; ++side) {
				fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.rTex[side].obj, 0);
				fbo.check();
				quadObj.draw({
					shader : burgersComputeFluxSlopeShader[side],
					uniforms : {
						dpos : [1/this.nx, 1/this.nx]
					},
					texs : [
						this.qTex, 
						this.solidTex,
						this.uiTex
					]
				});
				fbo.unbind();
			}

			//construct flux:
			for (var side = 0; side < 2; ++side) {
				fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.fluxTex[side].obj, 0);
				fbo.check();
				quadObj.draw({
					shader : burgersComputeFluxShader[this.fluxMethod][side],
					uniforms : {
						dpos : [1/this.nx, 1/this.nx],
						dt_dx : dt / dxi[side]
					},
					uniforms : {
						dpos : [1/this.nx, 1/this.nx],
					},
					texs : [
						this.qTex,
						this.solidTex,
						this.uiTex,
						this.rTex[side]
					]
				});
				fbo.unbind();
			}

			//update state
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.scratchTex.obj, 0);
			fbo.check();
			quadObj.draw({
				shader : burgersUpdateStateShader,
				uniforms : {
					dpos : [1/this.nx, 1/this.nx],
					dt_dx : [
						dt / dx,
						dt / dy
					]
				},
				texs : [
					this.qTex, 
					this.solidTex,
					this.fluxTex[0], 
					this.fluxTex[1]
				]
			});
			fbo.unbind();
			this.swapTexs('scratchTex', 'qTex');
		
			//boundary again
			this.boundary();

			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.pressureTex.obj, 0);
			fbo.check();
			//compute pressure
			quadObj.draw({
				shader : burgersComputePressureShader,
				uniforms : {
					dpos : [1/this.nx, 1/this.nx],
					rangeMin : [xmin, ymin],
					rangeMax : [xmax, ymax],
					externalForce : [externalForceX, externalForceY],
					gamma : this.gamma
				},
				texs : [this.qTex, this.solidTex]
			});
			fbo.unbind();

			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.scratchTex.obj, 0);
			fbo.check();
			//apply momentum diffusion
			quadObj.draw({
				shader : burgersApplyPressureToMomentumShader,
				uniforms : {
					dpos : [1/this.nx, 1/this.nx],
					dt : dt,
					dt_dx : [dt / dx, dt / dy],
					externalForce : [externalForceX, externalForceY]
				},
				texs : [this.qTex, this.solidTex, this.pressureTex]
			});
			fbo.unbind();
			this.swapTexs('scratchTex', 'qTex');

			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.scratchTex.obj, 0);
			fbo.check();
			//apply work diffusion
			quadObj.draw({
				shader : burgersApplyPressureToWorkShader,
				uniforms : {
					dpos : [1/this.nx, 1/this.nx],
					dt : dt,
					dt_dx : [dt / dx, dt / dy],
					externalForce : [externalForceX, externalForceY]
				},
				texs : [this.qTex, this.solidTex, this.pressureTex]
			});
			fbo.unbind();
			this.swapTexs('scratchTex', 'qTex');	
		
		}
	},
	'Riemann / Roe' : {
		initStep : function() {
			for (var side = 0; side < 2; ++side) {
				fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.roeValueTex[side].obj, 0);
				fbo.check();
				quadObj.draw({
					shader : roeComputeRoeValueShader[side],
					texs : [this.qTex]
				});
				fbo.unbind();
			
				fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.eigenvalueTex[side].obj, 0);
				fbo.check();
				quadObj.draw({
					shader : roeComputeEigenvalueShader[side],
					texs : [this.qTex, this.roeValueTex[side]]
				});
				fbo.unbind();

				for (var column = 0; column < 4; ++column) {
					fbo.bind();
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.eigenvectorColumnTex[side][column].obj, 0);
					fbo.check();
					quadObj.draw({
						shader : roeComputeEigenvectorColumnShader[side][column],
						texs : [this.qTex, this.roeValueTex[side]]
					});
					fbo.unbind();
				}

				for (var column = 0; column < 4; ++column) {
					fbo.bind();
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.eigenvectorInverseColumnTex[side][column].obj, 0);
					fbo.check();
					quadObj.draw({
						shader : roeComputeEigenvectorInverseColumnShader[side][column],
						texs : [this.qTex, this.roeValueTex[side]]
					});
					fbo.unbind();
				}
			}
		},
		calcCFLTimestep : function() {
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.cflReduceTex.obj, 0);
			fbo.check();
			quadObj.draw({
				shader : roeComputeCFLShader,
				texs : [this.eigenvalueTex[0], this.eigenvalueTex[1]]
			});
			fbo.unbind();
			return this.reduceToDetermineCFL();
		},
		advect : function(dt) {
			var dx = (xmax - xmin) / this.nx;
			var dy = (ymax - ymin) / this.nx;
			var dxi = [dx, dy];
			
			for (var side = 0; side < 2; ++side) {
				fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.dqTildeTex[side].obj, 0);
				fbo.check();
				quadObj.draw({
					shader : roeComputeDeltaQTildeShader[side],
					texs : [
						this.qTex, 
						this.eigenvectorInverseColumnTex[side][0],
						this.eigenvectorInverseColumnTex[side][1],
						this.eigenvectorInverseColumnTex[side][2],
						this.eigenvectorInverseColumnTex[side][3]
					]
				});
				fbo.unbind();
			}

			for (var side = 0; side < 2; ++side) {
				fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.rTildeTex[side].obj, 0);
				fbo.check();
				quadObj.draw({
					shader : roeComputeFluxSlopeShader[side],
					texs : [
						this.dqTildeTex[side],
						this.eigenvalueTex[side]
					]
				});
				fbo.unbind();
			}
	
			for (var side = 0; side < 2; ++side) {
				fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.fluxTex[side].obj, 0);
				fbo.check();
				quadObj.draw({
					shader : roeComputeFluxShader[this.fluxMethod][side],
					uniforms : {
						dt_dx : dt / dxi[side]
					},
					texs : [
						this.qTex,
						this.dqTildeTex[side],
						this.rTildeTex[side],
						this.eigenvalueTex[side],
						this.eigenvectorInverseColumnTex[side][0],
						this.eigenvectorInverseColumnTex[side][1],
						this.eigenvectorInverseColumnTex[side][2],
						this.eigenvectorInverseColumnTex[side][3],
						this.eigenvectorColumnTex[side][0],
						this.eigenvectorColumnTex[side][1],
						this.eigenvectorColumnTex[side][2],
						this.eigenvectorColumnTex[side][3]
					]
				});
				fbo.unbind();
			}

			//update state
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.scratchTex.obj, 0);
			fbo.check();
			quadObj.draw({
				shader : roeUpdateStateShader,
				uniforms : {
					dt_dx : [
						dt / dx,
						dt / dy
					]
				},
				texs : [
					this.qTex, 
					this.solidTex,
					this.fluxTex[0], 
					this.fluxTex[1]
				]
			});
			fbo.unbind();
			this.swapTexs('scratchTex', 'qTex');
		}
	}
};

/*
EulerEquationBurgersSolver
	EulerEquationBurgersExplicit	<- getting carbuncle inaccuracies
	EulerEquationBurgersBackwardEulerGaussSeidel	<- not on the GPU yet
GodunovSolver	<- not stable yet
	EulerEquationGodunovSolver
		EulerEquationGodunovExplicit
		EulerEquationRoeExplicit
*/

var eulerEquationSimulation = {
	methods : {
	/*
		'Burgers / Explicit' : EulerEquationBurgersExplicit.prototype,
		'Burgers / Backward Euler via Gauss Seidel' : EulerEquationBurgersBackwardEulerGaussSeidel.prototype,
		'Godunov / Explicit' : EulerEquationGodunovExplicit.prototype,
		'Roe / Explicit' : EulerEquationRoeExplicit.prototype
	*/
	},
	initialConditions : {
		Sod : function() {
			var thiz = this;
			gl.viewport(0, 0, this.nx, this.nx);
			fbo.setColorAttachmentTex2D(0, this.qTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : resetSodShader,
						uniforms : {
							noiseAmplitude : useNoise ? noiseAmplitude : 0
						},
						texs : [thiz.noiseTex]
					});
				}
			});
			fbo.setColorAttachmentTex2D(0, this.solidTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : solidShader
					});
				}
			});
		},
		SodCylinder : function() {
			var thiz = this;
			gl.viewport(0, 0, this.nx, this.nx);
			fbo.setColorAttachmentTex2D(0, this.qTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : resetSodShader,
						uniforms : {
							noiseAmplitude : useNoise ? noiseAmplitude : 0
						},
						texs : [thiz.noiseTex]
					});
				}
			});
			fbo.setColorAttachmentTex2D(0, this.solidTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : resetSodCylinderSolidShader
					});
				}
			});
		},
		Wave : function() {
			var thiz = this;
			gl.viewport(0, 0, this.nx, this.nx);
			fbo.setColorAttachmentTex2D(0, this.qTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : resetWaveShader,
						uniforms : {
							noiseAmplitude : useNoise ? noiseAmplitude : 0
						},
						texs : [thiz.noiseTex]
					});
				}
			});
			fbo.setColorAttachmentTex2D(0, this.solidTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : solidShader
					});
				}
			});
		},
		//http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
		KelvinHemholtz : function() {
			var thiz = this;
			gl.viewport(0, 0, this.nx, this.nx);
			fbo.setColorAttachmentTex2D(0, this.qTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : resetKelvinHemholtzShader,
						uniforms : {
							noiseAmplitude : useNoise ? noiseAmplitude : 0
						},
						texs : [thiz.noiseTex]
					});
				}
			});
			fbo.setColorAttachmentTex2D(0, this.solidTex);
			fbo.draw({
				callback : function() {
					quadObj.draw({
						shader : solidShader
					});
				}
			});
		}
	}
};

var HydroState = makeClass({ 
	init : function(args) {
		var thiz = this;

		this.nx = args.size;
		this.cfl =.5;
		this.gamma = args.gamma;

		var dpos = [1/this.nx, 1/this.nx];
		
		//provide default dpos values to allShaders
		$.each(allShaders, function(i,shader) {
			shader
				.use()
				.setUniforms({
					dpos : dpos,
					gamma : thiz.gamma,
					rangeMin : [xmin, ymin],
					rangeMax : [xmax, ymax],
					externalForce : [externalForceX, externalForceY]
				})
				.useNone();
		});

		//http://lab.concord.org/experiments/webgl-gpgpu/webgl.html
		encodeTempTex = new glutil.Texture2D({
			internalFormat : gl.RGBA,
			format : gl.RGBA,
			type : gl.UNSIGNED_BYTE,
			width : this.nx,
			height : this.nx,
			minFilter : gl.NEAREST,
			magFilter : gl.NEAREST,
			wrap : {
				s : gl.REPEAT,
				t : gl.REPEAT
			}
		});

		this.noiseTex = new glutil.Texture2D({
			internalFormat : gl.RGBA,
			format : gl.RGBA,
			type : gl.FLOAT,
			width : this.nx,
			height : this.nx,
			minFilter : gl.NEAREST,
			magFilter : gl.NEAREST,
			wrap : {
				s : gl.REPEAT,
				t : gl.REPEAT
			},
			data : function(i,j) {
				return [
					Math.random(),
					Math.random(),
					Math.random(),
					Math.random()
				];
			}
		});

		this.allFloatTexs = [];


		//scratch target used for any double-buffered operations
		this.allFloatTexs.push(this.scratchTex = new FloatTexture2D(this.nx, this.nx));

		//I'm skipping on x_i,j
		//Instead just use the uniforms xmin, xmax, ymin, ymax

		//I'm skipping on x_{i-1/2,j-1/2} because 
		//(1) you can reconstruct it from x
		//(2) it is only used for calculating dx's
		
		//q_i,j,state: state vector, stored as q[state + 4 * (j + this.nx * i)]
		//q_i,j,0: density: rho
		//q_i,j,1: momentum: rho * u
		//q_i,j,2: momentum: rho * v
		//q_i,j,3: work: rho * e
		this.qTex = new FloatTexture2D(this.nx, this.nx);	//rho, rho * u, rho * v, rho * e
		this.allFloatTexs.push(this.qTex);
		
		this.allFloatTexs.push(this.solidTex = new FloatTexture2D(this.nx, this.nx));

		//p_i,j: pressure
		this.pressureTex = new FloatTexture2D(this.nx, this.nx);
		this.allFloatTexs.push(this.pressureTex);

		//TODO it is tempting to merge r, f, and ui into an edge structure
		//and associate them with the nodes on either side of them,
		//but then I would lose out on the 2nd-order contributions to the flux limiter.

		//f_{i-1/2},{j-1/2},side,state: cell flux
		this.fluxTex = [];
		for (var side = 0; side < 2; ++side) {
			this.fluxTex[side] = new FloatTexture2D(this.nx, this.nx);
			this.allFloatTexs.push(this.fluxTex[side]);
		}

		this.cflReduceTex = new FloatTexture2D(this.nx, this.nx);
		this.allFloatTexs.push(this.cflReduceTex);
		this.nextCFLReduceTex = new FloatTexture2D(this.nx, this.nx);
		this.allFloatTexs.push(this.nextCFLReduceTex);


		//used for Burgers
		
		
		//r_{i-1/2},{j-1/2},side,state	
		this.rTex = [];
		for (var side = 0; side < 2; ++side) {
			this.rTex[side] = new FloatTexture2D(this.nx, this.nx);
			this.allFloatTexs.push(this.rTex[side]);
		}
		
		//only used with Burger's eqn advection code
		//u_{i-1/2},{j-1/2},dim: interface velocity
		this.uiTex = new FloatTexture2D(this.nx, this.nx);
		this.allFloatTexs.push(this.uiTex);


		//used for Riemann
	

		//v_{i-1/2},{j-1/2},side = [velocity.x, velocity.y, hTotal, speedOfSound]
		this.roeValueTex = [];
		for (var side = 0; side < 2; ++side) {
			this.roeValueTex[side] = new FloatTexture2D(this.nx, this.nx);
			this.allFloatTexs.push(this.roeValueTex[side]);
		}

		//a_{i-1/2},{j-1/2},side,state,state
		this.eigenvalueTex  = [];
		for (var side = 0; side < 2; ++side) {
			this.eigenvalueTex[side] = new FloatTexture2D(this.nx, this.nx); 
			this.allFloatTexs.push(this.eigenvalueTex[side]);
		}
		this.eigenvectorColumnTex = [];
		for (var side = 0; side < 2; ++side) {
			this.eigenvectorColumnTex[side] = [];
			for (var state = 0; state < 4; ++state) {
				this.eigenvectorColumnTex[side][state] = new FloatTexture2D(this.nx, this.nx);
				this.allFloatTexs.push(this.eigenvectorColumnTex[side][state]);
			}
		}
		this.eigenvectorInverseColumnTex = [];
		for (var side = 0; side < 2; ++side) {
			this.eigenvectorInverseColumnTex[side] = [];
			for (var state = 0; state < 4; ++state) {
				this.allFloatTexs.push(this.eigenvectorInverseColumnTex[side][state] = new FloatTexture2D(this.nx, this.nx));
			}
		}
		
		//qiTilde_{i-1/2},{j-1/2},side,state	
		this.dqTildeTex = [];
		for (var side = 0; side < 2; ++side) {
			this.allFloatTexs.push(this.dqTildeTex[side] = new FloatTexture2D(this.nx, this.nx));
		}

		//rTilde_{i-1/2},{j-1/2},side,state
		this.rTildeTex = []; 
		for (var side = 0; side < 2; ++side) {
			this.allFloatTexs.push(this.rTildeTex[side] = new FloatTexture2D(this.nx, this.nx)); 
		}

		//solver configuration
		this.boundaryMethod = 'periodic';
		this.fluxMethod = 'superbee';
		this.advectMethod = 'Riemann / Roe';
		this.drawToScreenMethod = 'Density';
	
		//initialize all textures to zero by default
		//...except the noise tex, which isn't added to the 'allFloatTexs' list
		gl.viewport(0, 0, this.nx, this.nx);
		fbo.bind();
		$.each(this.allFloatTexs, function(i,tex) {
			fbo.setColorAttachmentTex2D(0, tex);
			fbo.check();
			quadObj.draw({
				shader : solidShader
			});
		});
		fbo.unbind();

		//initial conditions
		this.resetSod();
		//eulerEquations.initialConditions.Sod.call(this);
	},
	//begin delete
	resetSod : function() {
		var thiz = this;
		gl.viewport(0, 0, this.nx, this.nx);
		fbo.setColorAttachmentTex2D(0, this.qTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : resetSodShader,
					uniforms : {
						noiseAmplitude : useNoise ? noiseAmplitude : 0
					},
					texs : [thiz.noiseTex]
				});
			}
		});
		fbo.setColorAttachmentTex2D(0, this.solidTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : solidShader
				});
			}
		});
	},
	resetSodCylinder : function() {
		var thiz = this;
		gl.viewport(0, 0, this.nx, this.nx);
		fbo.setColorAttachmentTex2D(0, this.qTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : resetSodShader,
					uniforms : {
						noiseAmplitude : useNoise ? noiseAmplitude : 0
					},
					texs : [thiz.noiseTex]
				});
			}
		});
		fbo.setColorAttachmentTex2D(0, this.solidTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : resetSodCylinderSolidShader
				});
			}
		});
	},
	resetWave : function() {
		var thiz = this;
		gl.viewport(0, 0, this.nx, this.nx);
		fbo.setColorAttachmentTex2D(0, this.qTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : resetWaveShader,
					uniforms : {
						noiseAmplitude : useNoise ? noiseAmplitude : 0
					},
					texs : [thiz.noiseTex]
				});
			}
		});
		fbo.setColorAttachmentTex2D(0, this.solidTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : solidShader
				});
			}
		});
	},
	//http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
	resetKelvinHemholtz : function() {
		var thiz = this;
		gl.viewport(0, 0, this.nx, this.nx);
		fbo.setColorAttachmentTex2D(0, this.qTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : resetKelvinHemholtzShader,
					uniforms : {
						noiseAmplitude : useNoise ? noiseAmplitude : 0
					},
					texs : [thiz.noiseTex]
				});
			}
		});
		fbo.setColorAttachmentTex2D(0, this.solidTex);
		fbo.draw({
			callback : function() {
				quadObj.draw({
					shader : solidShader
				});
			}
		});
	},
	//end delete
	boundary : function() {
		boundaryMethods[this.boundaryMethod].call(this);
	},
	step : function(dt) {
		var dx = (xmax - xmin) / this.nx;
		var dy = (ymax - ymin) / this.nx;
		
		//apply boundary conditions
		this.boundary();

		//solve
		advectMethods[this.advectMethod].advect.call(this, dt);
	},
	update : function() {
		gl.viewport(0, 0, this.nx, this.nx);
		//fbo.bind();
			
		//do any pre-calcCFLTimestep preparation (Roe computes eigenvalues here)
		advectMethods[this.advectMethod].initStep.call(this);

		//get timestep
		var dt;
		if (useCFL) {
			dt = advectMethods[this.advectMethod].calcCFLTimestep.call(this);
		} else {
			dt = fixedDT;
		}
window.lastDT = dt;
		//do the update
		this.step(dt);
			
		//fbo.unbind();
	},

	swapTexs : function(texFieldA, texFieldB) {
		//swap
		var tmp = this[texFieldA];
		this[texFieldA] = this[texFieldB];
		this[texFieldB] = tmp; 
	},

	//reduce to determine CFL
	reduceToDetermineCFL : function() {
		var size = this.nx;
		while (size > 1) {
			//console.log(getFloatTexData({srcTex:this.cflReduceTex})); fbo.bind();
			//console.log('reducing...')
			
			size /= 2;
			if (size !== Math.floor(size)) throw 'got a npo2 size '+this.nx;
			gl.viewport(0, 0, size, size);
			
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.scratchTex.obj, 0);
			fbo.check();
			gl.clear(gl.COLOR_BUFFER_BIT);
			quadObj.draw({
				shader : minReduceShader,
				uniforms : {
					texsize : [this.nx, this.nx], 
					viewsize : [size, size]
				},
				texs : [this.cflReduceTex]
			});
			fbo.unbind();

			this.swapTexs('scratchTex', 'cflReduceTex');
		}
		
		//now that the viewport is 1x1, run the encode shader on it
		fbo.bind();
		gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, encodeTempTex.obj, 0);
		fbo.check();
		gl.viewport(0, 0, encodeTempTex.width, encodeTempTex.height);
		quadObj.draw({
			shader : encodeShader[0],
			texs : [this.cflReduceTex]
		});

		var cflUint8Result = new Uint8Array(4);
		gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, cflUint8Result);
		
		fbo.unbind();
		
		var cflFloat32Result = new Float32Array(cflUint8Result.buffer);
		gl.viewport(0, 0, this.nx, this.nx);
		var result = cflFloat32Result[0] * this.cfl;
		return result;
	},

	/*
	args:
		vtxs
		texCoords (optional)
		everything else is forwarded to SceneObject.draw
	*/
	drawLine : function(args) {
		for (var i = 0; i < args.vtxs.length; ++i) {
			args.vtxs[i] *= 2;
			args.vtxs[i] -= 1;
		}
		lineObj.attrs.vertex.data.set(args.vtxs);
		lineObj.attrs.vertex.updateData();
		if (args.texCoords !== undefined) lineObj.attrs.texCoord.data.set(args.texCoords);
		lineObj.attrs.texCoord.updateData();
		lineObj.draw(args);
	}
});


var Hydro = makeClass({
	init : function() {
		var size = Number($.url().param('size'));
		if (size === undefined || size !== size) size = 512;
		var gamma = Number($.url().param('gamma'));
		if (gamma === undefined || gamma !== gamma) gamma = 7/5;
		
		this.state = new HydroState({
			size : size,
			gamma : gamma
		});
	},
	update : function() {
		
		//todo adm or something
		//update a copy of the grid and its once-refined
		//...and a once-unrefined ... over mergeable cells only?
		//then test for errors and split when needed
		this.state.update();

		//TODO reduce shader to determine min and max of whatever rendered value we will be using
		//until then, fixed only!
	}
});

var hydro;

function update() {
	//iterate
	hydro.update();

	/*
	//extract eigenvector and inverse textures
	if (hydro.state.advectMethod == 'Riemann / Roe') {
		var nx = hydro.state.nx;
		var maxmaxerr = 0;
		for (var side = 0; side < 2; ++side) {
			var eigenvectors = [];
			var eigenvectorInverses = [];
			for (var row = 0; row < 4; ++row) {
				eigenvectors[row] = [];
				eigenvectorInverses[row] = [];
				for (var column = 0; column < 4; ++column) {
					eigenvectors[row][column] = getFloatTexData({srcTex:hydro.state.eigenvectorColumnTex[side][column], channel:row});
					eigenvectorInverses[row][column] = getFloatTexData({srcTex:hydro.state.eigenvectorInverseColumnTex[side][column], channel:row});
				}
			}
		
			//now test them all
			for (var i = 0; i < 4; ++i) {
				for (var j = 0; j < 4; ++j) {
					var maxErr = 0;
					for (var e = 0; e < nx * nx; ++e) {
						var s = 0;
						for (var k = 0; k < 4; ++k) {
							s += eigenvectors[i][k][e] * eigenvectorInverses[k][j][e];
						}
						if (i == j) s--;
						maxErr = Math.max(maxErr, Math.abs(s));
					}
					maxmaxerr = Math.max(maxmaxerr, maxErr);
					console.log('side',side,'coord',i,j,'maxerr',maxErr);
				}
			}
		}
		console.log('greatest max err',maxmaxerr);
		console.log('done!');
		hydro.state.advectMethod = 'Burgers';
	//looks like we're close to the identity matrix, but we're getting some error, and it's building in the system
	}
	*/

	//reset viewport
	gl.viewport(0, 0, glutil.canvas.width, glutil.canvas.height);
	
	//draw
	glutil.draw();

	hydro.state.qTex.bind(0);
	if (gl.getExtension('OES_texture_float_linear')) gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
	hydro.state.pressureTex.bind(1);
	if (gl.getExtension('OES_texture_float_linear')) gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
	currentColorScheme.bind(2);

	glutil.unitQuad.draw({
		shader : drawToScreenShader[hydro.state.drawToScreenMethod],
		uniforms : {
			lastMin : hydro.lastDataMin,
			lastMax : hydro.lastDataMax
		}
	});

	gl.bindTexture(gl.TEXTURE_2D, null);

	gl.activeTexture(gl.TEXTURE1);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
	gl.bindTexture(gl.TEXTURE_2D, null);

	gl.activeTexture(gl.TEXTURE0);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
	gl.bindTexture(gl.TEXTURE_2D, null);

	requestAnimFrame(update);
}

function onresize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	$('#content').height(window.innerHeight - 50);
	glutil.resize();
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

/*
args:
	fbo : fbo (default window.fbo)
	srcTex : srcTex
	destTex : destTex - should be RGBA / UNSIGNED_BYTE (default window.encodeTempTex)
	channel : channel (default 0)
*/
function getFloatTexData(args) {//fbo, srcTex, destTex, channel) {
	var fbo = args.fbo === undefined ? window.fbo : args.fbo;
	var srcTex = args.srcTex;
	var destTex = args.destTex === undefined ? window.encodeTempTex : args.destTex;
	var channel = args.channel === undefined ? 0 : args.channel;

	var destUint8Array = new Uint8Array(destTex.width * destTex.height * 4);
	fbo.setColorAttachmentTex2D(0, destTex);
	fbo.draw({
		callback : function() {
			gl.viewport(0, 0, destTex.width, destTex.height);
			quadObj.draw({
				shader : encodeShader[channel],
				texs : [srcTex]
			});
			gl.readPixels(0, 0, destTex.width, destTex.height, destTex.format, destTex.type, destUint8Array);
		}
	});
	var destFloat32Array = new Float32Array(destUint8Array.buffer);
window.destUint8Array = destUint8Array;
window.destFloat32Array = destFloat32Array;
	return destFloat32Array;
}


var sceneObjects = [];

var currentColorScheme;
var colorSchemes = {};

$(document).ready(function(){
	
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()
	
	//try {
		glutil = new GLUtil({canvas:canvas});
		gl = glutil.context;
	/*} catch (e) {
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}*/

	FloatTexture2D = makeClass({
		super : glutil.Texture2D,
		init : function(width, height) {
			var args = {};
			args.width = width;
			args.height = height;
			args.internalFormat = gl.RGBA;
			args.format = gl.RGBA;
			args.type = gl.FLOAT;
			args.minFilter = gl.NEAREST;
			args.magFilter = gl.NEAREST;
			args.wrap = {
				s : gl.REPEAT,
				t : gl.REPEAT
			};
			FloatTexture2D.super.call(this, args);
		}
	});

	lineObj = new glutil.SceneObject({
		scene : glutil.scene,
		mode : gl.LINES,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				dim : 2,
				data : [0, 0, 1, 1],
				usage : gl.DYNAMIC_DRAW
			}),
			texCoord : new glutil.ArrayBuffer({
				dim : 2,
				data : [0, 0, 1, 1],
				usage : gl.DYNAMIC_DRAW
			})
		},
		parent : null,
		static : true
	});

	quadObj = new glutil.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				dim : 2,
				data : [-1,-1, 1,-1, -1,1, 1,1]
			}),
			texCoord : new glutil.ArrayBuffer({
				dim : 2,
				data : [0,0, 1,0, 0,1, 1,1]
			})
		},
		parent : null,
		static : true
	});

	//init shaders before init hydro (so it can call the resetSod or whatever)
	
		resetSodShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	vec2 gridPos = rangeMin + pos * (rangeMax - rangeMin);
	//I would use dpos functions, but it's only for initialization...
	float rho;
	if (gridPos.x < (.7 * rangeMin.x + .3 * rangeMax.x) 
		&& gridPos.y < (.7 * rangeMin.y + .3 * rangeMax.y))
	{
		rho = 1.;
	} else {
		rho = .1;
	}
	vec2 vel = vec2(0., 0.);
	vel.xy += (texture2D(randomTex, pos).xy - .5) * 2. * noiseAmplitude;
	float energyKinetic = .5 * dot(vel, vel);
	float energyPotential = dot(gridPos - rangeMin, externalForce);
	float energyThermal = 1.;
	float energyTotal = energyKinetic + energyThermal + energyPotential;
	gl_FragColor = vec4(rho, rho * vel.x, rho * vel.y, rho * energyTotal);
}		
*/}),
			uniforms : {
				rangeMin : ['vec2', [xmin, ymin]],
				rangeMax : ['vec2', [xmax, ymax]],
				externalForce : 'vec2',
				noiseAmplitude : 'float'
			},
			texs : ['randomTex']
		});

		resetSodCylinderSolidShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	vec2 gridPos = rangeMin + pos * (rangeMax - rangeMin);
	vec2 center = .35 * rangeMin + .65 * rangeMax;
	vec2 delta = gridPos - center;
	float distSq = dot(delta, delta);
	if (distSq < .1 * .1) {
		gl_FragColor = vec4(1., 1., 1., 1.);
	} else {
		gl_FragColor = vec4(0., 0., 0., 0.);
	}
}		
*/}),
			uniforms : {
				rangeMin : ['vec2', [xmin, ymin]],
				rangeMax : ['vec2', [xmax, ymax]],
			},
			texs : ['randomTex']
		});

		resetWaveShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	vec2 gridPos = rangeMin + pos * (rangeMax - rangeMin);
	float dg = .2 * (rangeMax.x - rangeMin.x);
	vec2 rangeMid = .5 * (rangeMax + rangeMin);
	vec2 gridPosFromCenter = gridPos - .5 * rangeMid;
	//Note I turned this down from 3. to 2. because GPU is currently using fixed dpos sizes
	float rho = 2. * exp(-dot(gridPosFromCenter, gridPosFromCenter) / (dg * dg)) + .1;
	vec2 vel = vec2(0., 0.);
	vel.xy += (texture2D(randomTex, pos).xy - .5) * 2. * noiseAmplitude;
	float energyKinetic = .5 * dot(vel, vel);
	float energyPotential = dot(gridPos - rangeMin, externalForce);
	float energyThermal = 1.;
	float energyTotal = energyKinetic + energyThermal + energyPotential;
	gl_FragColor = vec4(rho, rho * vel.xy, rho * energyTotal);
}		
*/}),
			uniforms : {
				rangeMin : ['vec2', [xmin, ymin]],
				rangeMax : ['vec2', [xmax, ymax]],
				externalForce : 'vec2',
				noiseAmplitude : 'float'
			},
			texs : ['randomTex']
		});

		resetKelvinHemholtzShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	vec2 gridPos = rangeMin + pos * (rangeMax - rangeMin);
	float rho = 1.;
	vec2 vel = vec2(0., 0.);
	if (gridPos.y > (.75 * rangeMin.y + .25 * rangeMax.y) && 
		gridPos.y < (.25 * rangeMin.y + .75 * rangeMax.y))
	{
		rho = 2.;
		vel.x = .5;
	} else {
		rho = 1.;
		vel.x = -.5;
	}
	vel.xy += (texture2D(randomTex, pos).xy - .5) * 2. * noiseAmplitude;
	//P = (gamma - 1) rho (eTotal - eKinetic)
	//eTotal = P / ((gamma - 1) rho) + eKinetic
	float pressure = 2.5;
	float energyKinetic = .5 * dot(vel, vel);
	float energyPotential = dot(gridPos - rangeMin, externalForce);
	float energyTotal = pressure / ((gamma - 1.) * rho) + energyKinetic + energyPotential; 
	gl_FragColor = vec4(rho, rho * vel.xy, rho * energyTotal);
}		
*/}),
			uniforms : {
				rangeMin : ['vec2', [xmin, ymin]],
				rangeMax : ['vec2', [xmax, ymax]],
				externalForce : 'vec2',
				noiseAmplitude : 'float',
				gamma : 'float'
			},
			texs : ['randomTex']
			//TODO make it periodic on the left/right borders and reflecting on the top/bottom borders	
		});

		burgersComputeCFLShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	float solid = texture2D(solidTex, pos).x;
	if (solid > .5) {
		//min? assign the max.
		gl_FragColor = vec4(1.);
	} else {
		vec2 gridPos = rangeMin + pos * (rangeMax - rangeMin);
		vec4 q = texture2D(qTex, pos);
		float rho = q.x;
		vec2 vel = q.yz / rho;
		float energyTotal = q.w / rho;
		float energyKinetic = .5 * dot(vel, vel);
		float energyPotential = dot(gridPos - rangeMin, externalForce);
		float energyThermal = energyTotal - energyKinetic - energyPotential;
		float speedOfSound = sqrt(gamma * (gamma - 1.) * energyThermal);
		vec2 dx = (rangeMax - rangeMin) * dpos;
		float cflx = dx.x / (speedOfSound + abs(vel.x));
		float cfly = dx.y / (speedOfSound + abs(vel.y));
		float cfl = min(cflx, cfly);
		gl_FragColor = vec4(cfl, 0., 0., 0.);
	}
}
*/}),
			uniforms : {
				rangeMin : ['vec2', [xmin, ymin]],
				rangeMax : ['vec2', [xmax, ymax]],
				dpos : 'vec2',
				externalForce : 'vec2',
				gamma : 'float'
			},
			texs : ['qTex', 'solidTex']
		});

		burgersComputeInterfaceVelocityShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	vec2 dposx = vec2(dpos.x, 0.);
	vec2 dposy = vec2(0., dpos.y);
	vec4 qR = texture2D(qTex, pos);
	vec4 qXL = texture2D(qTex, pos - dposx);
	vec4 qYL = texture2D(qTex, pos - dposy);
	float solidR = texture2D(solidTex, pos).x;
	float solidXL = texture2D(solidTex, pos - dposx).x;
	float solidYL = texture2D(solidTex, pos - dposy).x;
	float uR = qR.y / qR.x;
	float vR = qR.z / qR.x;
	float uL = qXL.y / qXL.x;
	float vL = qYL.z / qYL.x;
	if (solidXL > .5 && solidR <= .5) uL = -uR;
	if (solidR > .5 && solidXL <= .5) uR = -uL;
	if (solidYL > .5 && solidR <= .5) vL = -vR;
	if (solidR > .5 && solidYL <= .5) vR = -vL;
	gl_FragColor = vec4(
		.5 * (uR + uL),
		.5 * (vR + vL),
		0., 0.);
}		
*/}),
			uniforms : {
				dpos : 'vec2'
			},
			texs : ['qTex', 'solidTex']
		});

		$.each(coordNames, function(i, coordName) {
			burgersComputeFluxSlopeShader[i] = new glutil.KernelShader({
				code : mlstr(function(){/*
void main() {
	vec2 sidestep = vec2(0., 0.);
	sidestep[$side] = dpos[$side];

	vec2 posL2 = pos - 2.*sidestep;
	vec2 posL1 = pos - sidestep;
	vec2 posR1 = pos;
	vec2 posR2 = pos + sidestep;

	float solidL2 = texture2D(solidTex, posL2).x;
	float solidL1 = texture2D(solidTex, posL1).x;
	float solidR1 = texture2D(solidTex, posR1).x;
	float solidR2 = texture2D(solidTex, posR2).x;

	vec4 qL2 = texture2D(qTex, posL2);
	vec4 qL1 = texture2D(qTex, posL1);
	vec4 qR1 = texture2D(qTex, posR1);
	vec4 qR2 = texture2D(qTex, posR2);

	if (solidL2 > .5) {
		qL2 = qL1;
		qL2[$side+1] = -qL2[$side+1];
	}
	if (solidR2 > .5) {
		qR2 = qR1;
		qR2[$side+1] = -qR2[$side+1];
	}
	if (solidL1 > .5) {
		qL1 = -qR1;
		qL1[$side+1] = -qL1[$side+1];
	}
	if (solidR1 > .5) {
		qR1 = -qL1;
		qR1[$side+1] = -qR1[$side+1];
	}
	
	//dq = q_i,j - q_{{i,j}-dirs[side]}
	vec4 dq = qR1 - qL1; 
	
	vec4 s = sign(dq);
	s *= s;
	
	float ui = texture2D(uiTex, pos)[$side];
	float uigz = step(0., ui); 
	gl_FragColor = s * mix(qR2 - qR1, qL1 - qL2, uigz) / dq;
}
*/}).replace(/\$side/g, i),
				uniforms : {
					dpos : 'vec2'
				},
				texs : ['qTex', 'solidTex', 'uiTex']
			});
		});

		$.each(fluxMethods, function(methodName,fluxMethodCode) {
			burgersComputeFluxShader[methodName] = [];
			$.each(coordNames, function(i, coordName) {
				burgersComputeFluxShader[methodName][i] = new glutil.KernelShader({
					code : mlstr(function(){/*
vec4 fluxMethod(vec4 r) {
	$fluxMethodCode
}			
*/}).replace(/\$fluxMethodCode/g, fluxMethodCode)
+ mlstr(function(){/*
void main() {
	vec2 sidestep = vec2(0., 0.);
	sidestep[$side] = dpos[$side];
	
	vec2 posL = pos - sidestep;
	vec2 posR = pos;

	float ui = texture2D(uiTex, pos)[$side];

	vec4 qL = texture2D(qTex, posL);
	vec4 qR = texture2D(qTex, posR);

	float solidL = texture2D(solidTex, posL).x;
	float solidR = texture2D(solidTex, posR).x;

	if (solidL > .5 && solidR <= .5) {
		qL = qR;
		qL[$side+1] = -qL[$side+1];
	}
	if (solidR > .5 && solidL <= .5) {
		qR = qL;
		qR[$side+1] = -qR[$side+1];
	}

	if (ui >= 0.) {
		gl_FragColor = ui * qL;
	} else {
		gl_FragColor = ui * qR;
	}
	
	vec4 r = texture2D(rTex, pos);
	vec4 phi = fluxMethod(r);
	vec4 delta = phi * (qR - qL);
	
	gl_FragColor += delta * .5 * .5 * abs(ui) * (1. - abs(ui * dt_dx));
}			
*/}).replace(/\$side/g, i),
					uniforms : {
						dpos : 'vec2',
						dt_dx : 'float'
					},
					texs : ['qTex', 'solidTex', 'uiTex', 'rTex']
				});
			});
		});

		burgersUpdateStateShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	float solid = texture2D(solidTex, pos).x;
	if (solid > .5) {
		gl_FragColor = vec4(0., 0., 0., 0.);
	} else {
		vec4 q = texture2D(qTex, pos);
		vec4 fluxXL = texture2D(fluxXTex, pos);
		vec4 fluxXR = texture2D(fluxXTex, pos + vec2(dpos.x, 0.));
		vec4 fluxYL = texture2D(fluxYTex, pos);
		vec4 fluxYR = texture2D(fluxYTex, pos + vec2(0., dpos.y));
		gl_FragColor = q 
			- dt_dx.x * (fluxXR - fluxXL)
			- dt_dx.y * (fluxYR - fluxYL);
	}
}		
*/}),
			uniforms : {
				dpos : 'vec2',
				dt_dx : 'vec2'
			},
			texs : ['qTex', 'solidTex', 'fluxXTex', 'fluxYTex']
		});


		//used for Riemann

		roeComputeCFLShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	vec4 eigenvalueXN = texture2D(eigenvalueXTex, pos);
	vec4 eigenvalueXP = texture2D(eigenvalueXTex, pos + vec2(dpos.x, 0.));
	float maxLambdaXN = max(0., max(max(eigenvalueXN.x, eigenvalueXN.y), max(eigenvalueXN.z, eigenvalueXN.w)));
	float minLambdaXP = min(0., min(min(eigenvalueXP.x, eigenvalueXP.y), min(eigenvalueXP.z, eigenvalueXP.w)));
	
	vec4 eigenvalueYN = texture2D(eigenvalueYTex, pos);
	vec4 eigenvalueYP = texture2D(eigenvalueYTex, pos + vec2(0., dpos.y));
	float maxLambdaYN = max(0., max(max(eigenvalueYN.x, eigenvalueYN.y), max(eigenvalueYN.z, eigenvalueYN.w)));
	float minLambdaYP = min(0., min(min(eigenvalueYP.x, eigenvalueYP.y), min(eigenvalueYP.z, eigenvalueYP.w)));
	
	vec2 dxi = (rangeMax - rangeMin) * dpos;
	gl_FragColor = vec4(
		min( dxi.x / (maxLambdaXN - minLambdaXP), dxi.y / (maxLambdaYN - minLambdaYP)),
		0., 0., 0.);
}
*/}),
			uniforms : {
				dpos : 'vec2',
				rangeMin : ['vec2', [xmin, ymin]],
				rangeMax : ['vec2', [xmax, ymax]]
			},
			texs : ['eigenvalueXTex', 'eigenvalueYTex']
		});

		$.each(coordNames, function(i,coordName) {
			roeComputeRoeValueShader[i] = new glutil.KernelShader({
				code : mlstr(function(){/*
void main() {
	vec2 sidestep = vec2(0., 0.);
	sidestep[$side] = dpos[$side];
	
	vec2 posL = pos - sidestep;
	vec2 gridposL = rangeMin + posL * (rangeMax - rangeMin);
	vec4 qL = texture2D(qTex, posL);
	float densityL = qL.x;
	float invDensityL = 1. / densityL;
	vec2 velocityL = qL.yz * invDensityL; 
	float energyTotalL = qL.w * invDensityL;
	float energyKineticL = .5 * dot(velocityL, velocityL);
	float energyPotentialL = dot(gridposL - rangeMin, externalForce);
	float energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
	float pressureL = (gamma - 1.) * densityL * energyThermalL;
	float speedOfSoundL = sqrt(gamma * pressureL * invDensityL);
	float hTotalL = energyTotalL + pressureL * invDensityL;
	float roeWeightL = sqrt(densityL);
	
	vec2 posR = pos;
	vec2 gridposR = rangeMin + posR * (rangeMax - rangeMin);
	vec4 qR = texture2D(qTex, posR);
	float densityR = qR.x; 
	float invDensityR = 1. / densityR;
	vec2 velocityR = qR.yz * invDensityR;
	float energyTotalR = qR.w * invDensityR; 
	float energyKineticR = .5 * dot(velocityR, velocityR);
	float energyPotentialR = dot(gridposR - rangeMin, externalForce);
	float energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
	float pressureR = (gamma - 1.) * densityR * energyThermalR;
	float speedOfSoundR = sqrt(gamma * pressureR * invDensityR);
	float hTotalR = energyTotalR + pressureR * invDensityR;
	float roeWeightR = sqrt(densityR);

	float invDenom = 1. / (roeWeightL + roeWeightR);
	vec2 velocity = (roeWeightL * velocityL + roeWeightR * velocityR) * invDenom;
	float hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) * invDenom;
	
	float velocitySq = dot(velocity, velocity);
	float speedOfSound = sqrt((gamma - 1.) * (hTotal - .5 * velocitySq));
	
	gl_FragColor = vec4(
		velocity,
		hTotal,
		speedOfSound);
}
*/}).replace(/\$side/g, i),
				uniforms : {
					dpos : 'vec2',
					rangeMin : 'vec2',
					rangeMax : 'vec2',
					externalForce : 'vec2',
					gamma : 'float'
				},
				texs : ['qTex']
			});
		});

		$.each(coordNames, function(i,coordName) {	
			roeComputeEigenvalueShader[i] = new glutil.KernelShader({
				code : mlstr(function(){/*
void main() {
	vec4 roeValues = texture2D(roeValueTex, pos);
	vec2 velocity = roeValues.xy;
	float hTotal = roeValues.z; 
	float speedOfSound = roeValues.w;

	vec2 normal = vec2(0., 0.);
	normal[$side] = 1.;
	vec2 tangent = vec2(-normal.y, normal.x);
	
	float velocityN = dot(velocity, normal);
	float velocityT = dot(velocity, tangent);
	float velocitySq = dot(velocity, velocity);
	
	//eigenvalues: min, mid, max
	gl_FragColor = vec4(
		velocityN - speedOfSound,
		velocityN,
		velocityN,
		velocityN + speedOfSound);
}
*/}).replace(/\$side/g, i),
				uniforms : {
					gamma : 'float'
				},
				texs : ['qTex', 'roeValueTex']
			});
		});

		var roeComputeEigenvectorColumnCode = [
			mlstr(function(){/*
	//min eigenvector
	gl_FragColor = vec4(
		1.,
		velocity.x - speedOfSound * normal.x,
		velocity.y - speedOfSound * normal.y,
		hTotal - speedOfSound * velocityN);
*/}),
			mlstr(function(){/*
	//mid eigenvector (normal)
	gl_FragColor = vec4(
		1.,
		velocity.x,
		velocity.y,
		.5 * velocitySq);
*/}),
			mlstr(function(){/*
	//mid eigenvector (tangent)
	gl_FragColor = vec4(
		0.,
		tangent.x,
		tangent.y,
		velocityT);
*/}),
			mlstr(function(){/*
	//max eigenvector
	gl_FragColor = vec4(
		1.,
		velocity.x + speedOfSound * normal.x,
		velocity.y + speedOfSound * normal.y,
		hTotal + speedOfSound * velocityN);
*/})
		];

		$.each(coordNames, function(i,coordName) {	
			roeComputeEigenvectorColumnShader[i] = [];
			$.each(roeComputeEigenvectorColumnCode, function(j,code) {
				roeComputeEigenvectorColumnShader[i][j] = new glutil.KernelShader({
					code : mlstr(function(){/*
void main() {
	vec4 roeValues = texture2D(roeValueTex, pos);
	vec2 velocity = roeValues.xy;
	float hTotal = roeValues.z; 
	float speedOfSound = roeValues.w;

	vec2 normal = vec2(0., 0.);
	normal[$side] = 1.;
	vec2 tangent = vec2(-normal.y, normal.x);
	
	float velocityN = dot(velocity, normal);
	float velocityT = dot(velocity, tangent);
	float velocitySq = dot(velocity, velocity);
	
	
*/}).replace(/\$side/g, i) + code + '\n}',
					uniforms : {
						gamma : 'float'
					},
					texs : ['qTex', 'roeValueTex']
				});
			});
		});

		/*
		rows are min, mid normal, mid tangent, max
		but I'm going to write these out in columns for easier reconstruction of the original matrix 
		... at least I think it'll be easier
		*/
		var roeComputeEigenvectorInverseColumnCode = [
			mlstr(function(){/*
	float invDenom = .5 / (speedOfSound * speedOfSound);
	gl_FragColor = vec4(
		(.5 * (gamma - 1.) * velocitySq + speedOfSound * velocityN) * invDenom,	//ei_0,0
		1. - (gamma - 1.) * velocitySq * invDenom,	//ei_1,0
		-velocityT, //ei_2,0
		(.5 * (gamma - 1.) * velocitySq - speedOfSound * velocityN) * invDenom);	//ei_3,0
*/}),
			mlstr(function(){/*
	float invDenom = .5 / (speedOfSound * speedOfSound);
	gl_FragColor = vec4(
		-(normal.x * speedOfSound + (gamma - 1.) * velocity.x) * invDenom,	//ei_0,1
		(gamma - 1.) * velocity.x * 2. * invDenom,	//ei_1,1
		tangent.x,	//ei_2,1
		(normal.x * speedOfSound - (gamma - 1.) * velocity.x) * invDenom);	//ei_3,1
*/}),
			mlstr(function(){/*
	float invDenom = .5 / (speedOfSound * speedOfSound);
	gl_FragColor = vec4(
		-(normal.y * speedOfSound + (gamma - 1.) * velocity.y) * invDenom,	//ei_0,2
		(gamma - 1.) * velocity.y * 2. * invDenom,	//ei_1,2
		tangent.y,	//ei_2,2
		(normal.y * speedOfSound - (gamma - 1.) * velocity.y) * invDenom);	//ei_3,2
*/}),
			mlstr(function(){/*
	float invDenom = .5 / (speedOfSound * speedOfSound);
	gl_FragColor = vec4(
		(gamma - 1.) * invDenom,	//ei_0,3
		-(gamma - 1.) * 2. * invDenom,	//ei_1,3
		0.,	//ei_2,3
		(gamma - 1.) * invDenom);	//ei_3,3
*/})
		];

		$.each(coordNames, function(i,coordName) {	
			roeComputeEigenvectorInverseColumnShader[i] = [];
			$.each(roeComputeEigenvectorInverseColumnCode, function(j,code) {
				roeComputeEigenvectorInverseColumnShader[i][j] = new glutil.KernelShader({
					code : mlstr(function(){/*
void main() {
	vec4 roeValues = texture2D(roeValueTex, pos);
	vec2 velocity = roeValues.xy;
	float hTotal = roeValues.z; 
	float speedOfSound = roeValues.w;

	vec2 normal = vec2(0., 0.);
	normal[$side] = 1.;
	vec2 tangent = vec2(-normal.y, normal.x);
	
	float velocityN = dot(velocity, normal);
	float velocityT = dot(velocity, tangent);
	float velocitySq = dot(velocity, velocity);
	
*/}).replace(/\$side/g, i) + code + '\n}',
					uniforms : {
						gamma : 'float'
					},
					texs : ['qTex', 'roeValueTex']
				});
			});
		});

		$.each(coordNames, function(i, coordName) {
			roeComputeDeltaQTildeShader[i] = new glutil.KernelShader({
				code : mlstr(function(){/*
void main() {
	vec2 sidestep = vec2(0., 0.);
	sidestep[$side] = dpos[$side];

	vec4 eigenvectorInverseCol0 = texture2D(eigenvectorInverseCol0Tex, pos);
	vec4 eigenvectorInverseCol1 = texture2D(eigenvectorInverseCol1Tex, pos);
	vec4 eigenvectorInverseCol2 = texture2D(eigenvectorInverseCol2Tex, pos);
	vec4 eigenvectorInverseCol3 = texture2D(eigenvectorInverseCol3Tex, pos);
	mat4 eigenvectorInverse = mat4(
		eigenvectorInverseCol0,
		eigenvectorInverseCol1,
		eigenvectorInverseCol2,
		eigenvectorInverseCol3);

	vec4 qPrev = texture2D(qTex, pos - sidestep);
	vec4 q = texture2D(qTex, pos);
	vec4 dq = q - qPrev;

	gl_FragColor = eigenvectorInverse * dq;
}
*/}).replace(/\$side/g, i),
				uniforms : {
					dpos : 'vec2'
				},
				texs : [
					'qTex',
					'eigenvectorInverseCol0Tex',
					'eigenvectorInverseCol1Tex',
					'eigenvectorInverseCol2Tex',
					'eigenvectorInverseCol3Tex'
				]
			});
		});	

		$.each(coordNames, function(i,coordName) {
			roeComputeFluxSlopeShader[i] = new glutil.KernelShader({
				code : (mlstr(function(){/*
void main() {
	vec2 sidestep = vec2(0., 0.);
	sidestep[$side] = dpos[$side];
	vec4 dqTildePrev = texture2D(dqTildeTex, pos - sidestep); 
	vec4 dqTilde = texture2D(dqTildeTex, pos); 
	vec4 dqTildeNext = texture2D(dqTildeTex, pos + sidestep); 
	vec4 eigenvalues = texture2D(eigenvalueTex, pos);
*/}) + [0,1,2,3].map(function(j) {
	return mlstr(function(){/*
	if (abs(dqTilde[$j]) > 0.) {
		if (eigenvalues[$j] > 0.) {
			gl_FragColor[$j] = dqTildePrev[$j] / dqTilde[$j];
		} else {
			gl_FragColor[$j] = dqTildeNext[$j] / dqTilde[$j];
		}
	} else {
		gl_FragColor[$j] = 0.;
	}
*/}).replace(/\$j/g, j);
	}).join('')
+ '}').replace(/\$side/g, i),			
				uniforms : {
					dpos : 'vec2'
				},
				texs : ['dqTildeTex', 'eigenvalueTex']
			});
		});

		$.each(fluxMethods, function(methodName,fluxMethodCode) {
			roeComputeFluxShader[methodName] = [];
			$.each(coordNames, function(i,coordName) {
				roeComputeFluxShader[methodName][i] = new glutil.KernelShader({
					code : mlstr(function(){/*
vec4 fluxMethod(vec4 r) {
	$fluxMethodCode
}			
*/}).replace(/\$fluxMethodCode/g, fluxMethodCode)
+ mlstr(function(){/*
void main() {
	vec2 sidestep = vec2(0., 0.);
	sidestep[$side] = dpos[$side];

	vec4 q = texture2D(qTex, pos);
	vec4 qPrev = texture2D(qTex, pos - sidestep);

	vec4 eigenvalues = texture2D(eigenvalueTex, pos);

	vec4 eigenvectorCol0 = texture2D(eigenvectorCol0Tex, pos);
	vec4 eigenvectorCol1 = texture2D(eigenvectorCol1Tex, pos);
	vec4 eigenvectorCol2 = texture2D(eigenvectorCol2Tex, pos);
	vec4 eigenvectorCol3 = texture2D(eigenvectorCol3Tex, pos);
	mat4 eigenvectorMat = mat4(
		eigenvectorCol0, 
		eigenvectorCol1, 
		eigenvectorCol2, 
		eigenvectorCol3);
	
	vec4 eigenvectorInverseCol0 = texture2D(eigenvectorInverseCol0Tex, pos);
	vec4 eigenvectorInverseCol1 = texture2D(eigenvectorInverseCol1Tex, pos);
	vec4 eigenvectorInverseCol2 = texture2D(eigenvectorInverseCol2Tex, pos);
	vec4 eigenvectorInverseCol3 = texture2D(eigenvectorInverseCol3Tex, pos);
	mat4 eigenvectorInverseMat = mat4(
		eigenvectorInverseCol0, 
		eigenvectorInverseCol1, 
		eigenvectorInverseCol2, 
		eigenvectorInverseCol3);

	vec4 qAvg = (q + qPrev) * .5;
	vec4 fluxAvgTilde = eigenvectorInverseMat * qAvg;	///matrix multiply
	fluxAvgTilde *= eigenvalues;						//per-component multiply

	vec4 dqTilde = texture2D(dqTildeTex, pos);
	vec4 rTilde = texture2D(rTildeTex, pos);
	vec4 phi = fluxMethod(rTilde);

	vec4 theta = step(0., eigenvalues) * 2. - 1.;
	vec4 epsilon = eigenvalues * dt_dx;
	vec4 deltaFluxTilde = eigenvalues * dqTilde;
	vec4 fluxTilde = fluxAvgTilde - .5 * deltaFluxTilde * (theta + .5 * phi * (epsilon - theta));
	gl_FragColor = eigenvectorMat * fluxTilde;
}
*/}).replace(/\$side/g, i),
					uniforms : {
						dpos : 'vec2',
						dt_dx : 'float'
					},
					texs : [
						'qTex',
						'dqTildeTex',
						'rTildeTex',
						'eigenvalueTex',
						'eigenvectorInverseCol0Tex',
						'eigenvectorInverseCol1Tex',
						'eigenvectorInverseCol2Tex',
						'eigenvectorInverseCol3Tex',
						'eigenvectorCol0Tex',
						'eigenvectorCol1Tex',
						'eigenvectorCol2Tex',
						'eigenvectorCol3Tex',
					]
				});
			});
		});

		//matches burgersUpdateStateShader
		roeUpdateStateShader = burgersUpdateStateShader;


		//pressure shaders
		

		burgersComputePressureShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	float solid = texture2D(solidTex, pos).x;
	if (solid > .5) {
		gl_FragColor = vec4(0., 0., 0., 0.);
	} else {
		vec2 gridPos = rangeMin + pos * (rangeMax - rangeMin);
		vec4 q = texture2D(qTex, pos);
		float rho = q.x;
		vec2 vel = q.yz / rho;
		float energyTotal = q.w / rho;
		float energyKinetic = .5 * dot(vel, vel);
		float energyPotential = dot(gridPos - rangeMin, externalForce);
		float energyThermal = energyTotal - energyKinetic - energyPotential;
		gl_FragColor = vec4(0., 0., 0., 0.);
		gl_FragColor.x = (gamma - 1.) * rho * energyThermal;
	}
}
*/}),
			uniforms : {
				dpos : 'vec2',
				rangeMin : 'vec2',
				rangeMax : 'vec2',
				externalForce : 'vec2',
				gamma : 'float'
			},
			texs : ['qTex', 'solidTex']
		});

		burgersApplyPressureToMomentumShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	float solid = texture2D(solidTex, pos).x;
	if (solid > .5) {
		gl_FragColor = vec4(0., 0., 0., 0.);
	} else {	
		vec2 dposx = vec2(dpos.x, 0.);
		vec2 dposy = vec2(0., dpos.y);
		
		vec2 posXL = pos - dposx;
		vec2 posXR = pos + dposx;
		vec2 posYL = pos - dposy;
		vec2 posYR = pos + dposy;

		float solidXL = texture2D(solidTex, posXL).x;
		float solidXR = texture2D(solidTex, posXR).x;
		float solidYL = texture2D(solidTex, posYL).x;
		float solidYR = texture2D(solidTex, posYR).x;

		float pressureC = texture2D(pressureTex, pos).x;
		float pressureXL = texture2D(pressureTex, posXL).x;
		if (solidXL > .5) pressureXL = pressureC;
		float pressureXR = texture2D(pressureTex, posXR).x;
		if (solidXR > .5) pressureXR = pressureC;
		float pressureYL = texture2D(pressureTex, posYL).x;
		if (solidYL > .5) pressureYL = pressureC;
		float pressureYR = texture2D(pressureTex, posYR).x;
		if (solidYR > .5) pressureYR = pressureC; 
	
		vec2 dPressure = .5 * vec2(pressureXR - pressureXL, pressureYR - pressureYL);

		vec4 q = texture2D(qTex, pos);
		float rho = q.x;
		gl_FragColor = q;
		gl_FragColor.yz -= dt_dx * dPressure + dt * rho * externalForce;
	}
}
*/}),
			uniforms : {
				dpos : 'vec2',
				dt : 'float',
				dt_dx : 'vec2',
				externalForce : 'vec2',
			},
			texs : ['qTex', 'solidTex', 'pressureTex']
		});
		
		burgersApplyPressureToWorkShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	float solid = texture2D(solidTex, pos).x;
	if (solid > .5) {
		gl_FragColor = vec4(0., 0., 0., 0.);
	} else {
		vec2 dposx = vec2(dpos.x, 0.);
		vec2 dposy = vec2(0., dpos.y);
		
		vec2 posXR = pos + dposx;
		vec2 posXL = pos - dposx;
		vec2 posYR = pos + dposy;
		vec2 posYL = pos - dposy;

		float solidXL = texture2D(solidTex, posXL).x;
		float solidXR = texture2D(solidTex, posXR).x;
		float solidYL = texture2D(solidTex, posYL).x;
		float solidYR = texture2D(solidTex, posYR).x;

		float pressureC = texture2D(pressureTex, pos).x;
		float pressureXL = texture2D(pressureTex, posXL).x;
		if (solidXL > .5) pressureXL = pressureC;
		float pressureXR = texture2D(pressureTex, posXR).x;
		if (solidXR > .5) pressureXR = pressureC;
		float pressureYL = texture2D(pressureTex, posYL).x;
		if (solidYL > .5) pressureYL = pressureC;
		float pressureYR = texture2D(pressureTex, posYR).x;
		if (solidYR > .5) pressureYR = pressureC; 

		vec3 rho_C = texture2D(qTex, pos).xyz;
		vec2 rho_u_XR = texture2D(qTex, posXR).xy;
		vec2 rho_u_XL = texture2D(qTex, posXL).xy;
		vec2 rho_v_YR = texture2D(qTex, posYR).xz;
		vec2 rho_v_YL = texture2D(qTex, posYL).xz;

		vec2 uvC = rho_C.yz / rho_C.x;
		float uXR = rho_u_XR.y / rho_u_XR.x;
		if (solidXR > .5) uXR = -uvC.x;
		float uXL = rho_u_XL.y / rho_u_XL.x;
		if (solidXL > .5) uXL = -uvC.x;
		float vYR = rho_v_YR.y / rho_v_YR.x;
		if (solidYR > .5) vYR = -uvC.y;
		float vYL = rho_v_YL.y / rho_v_YL.x;
		if (solidYL > .5) vYL = -uvC.y;
		
		vec2 dPressureTimesVelocity = .5 * vec2(
			pressureXR * uXR - pressureXL * uXL, 
			pressureYR * vYR - pressureYL * vYL);

		vec4 q = texture2D(qTex, pos);
		gl_FragColor = q;
		gl_FragColor.w -= dot(dt_dx, dPressureTimesVelocity) + dt * dot(q.yz, externalForce);
	}
}
*/}),
			uniforms : {
				dpos : 'vec2',
				dt : 'float',
				dt_dx : 'vec2',
				externalForce : 'vec2'
			},
			texs : ['qTex', 'solidTex', 'pressureTex']
		});

		solidShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	gl_FragColor = color;
}		
*/}),
			uniforms : {
				color : ['vec4', [0,0,0,0]]
			}
		});

		copyShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	gl_FragColor = texture2D(srcTex, pos + offset);
}		
*/}),
			uniforms : {
				offset : ['vec2', [0,0]]
			},
			texs : ['srcTex']
		});

		minReduceShader = new glutil.KernelShader({
			code : mlstr(function(){/*
void main() {
	vec2 intPos = pos * viewsize - .5;
	
	float a = texture2D(srcTex, (intPos * 2. + .5) / texsize).x;
	float b = texture2D(srcTex, (intPos * 2. + vec2(1., 0.) + .5) / texsize).x;
	float c = texture2D(srcTex, (intPos * 2. + vec2(0., 1.) + .5) / texsize).x;
	float d = texture2D(srcTex, (intPos * 2. + vec2(1., 1.) + .5) / texsize).x;
	float e = min(a,b);
	float f = min(c,d);
	float g = min(e,f);
	gl_FragColor = vec4(g, 0., 0., 0.);
}
*/}),
			uniforms : {
				texsize : 'vec2',
				viewsize : 'vec2'
			},
			texs : ['srcTex']
		});

		fbo = new glutil.Framebuffer({
			width : this.nx,
			height : this.nx
		});
	
	var drawToScreenVertexShader = new glutil.VertexShader({
		code : glutil.vertexPrecision + mlstr(function(){/*
attribute vec2 vertex;
varying vec2 pos; 
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	pos = vertex;
	gl_Position = projMat * mvMat * vec4(vertex.xy, 0., 1.);
}
*/})
	});

	$.each(drawToScreenMethods, function(drawToScreenMethodName, drawToScreenMethodCode) {
		drawToScreenShader[drawToScreenMethodName] = new glutil.ShaderProgram({
			vertexShader : drawToScreenVertexShader,
			fragmentCode : mlstr(function(){/*
#extension GL_OES_standard_derivatives : enable
varying vec2 pos;
uniform vec2 dpos;
uniform vec2 rangeMin;
uniform vec2 rangeMax;
uniform float gamma;
uniform float lastMin, lastMax;
uniform sampler2D qTex;
uniform sampler2D pressureTex;
uniform sampler2D gradientTex;
*/}) + mlstr(function(){/*
float drawToScreenMethod() {
	$drawToScreenMethodCode
}			
*/}).replace(/\$drawToScreenMethodCode/g, drawToScreenMethodCode)
+ mlstr(function(){/*
void main() {
	float v = (drawToScreenMethod() - lastMin) / (lastMax - lastMin);
	gl_FragColor = texture2D(gradientTex, vec2(v, .5)); 
}	
*/}),
			fragmentPrecision : 'best',
			uniforms : {
				qTex : 0,
				pressureTex : 1,
				gradientTex : 2
			}
		});
	});


	//burgers
	allShaders.push(burgersComputeCFLShader);
	allShaders.push(burgersComputeInterfaceVelocityShader);
	allShaders = allShaders.concat(burgersComputeFluxSlopeShader);	
	$.each(burgersComputeFluxShader, function(k, fluxShaders) {
		allShaders = allShaders.concat(fluxShaders);
	});
	allShaders.push(burgersUpdateStateShader);
	allShaders.push(burgersComputePressureShader);
	allShaders.push(burgersApplyPressureToMomentumShader);
	allShaders.push(burgersApplyPressureToWorkShader);

	//riemann /roe
	allShaders.push(roeComputeCFLShader);
	allShaders = allShaders.concat(roeComputeRoeValueShader);
	allShaders = allShaders.concat(roeComputeEigenvalueShader);
	$.each(roeComputeEigenvectorColumnShader, function(k, evShaders) {
		allShaders = allShaders.concat(evShaders);
	});
	$.each(roeComputeEigenvectorInverseColumnShader, function(k, evInvShaders) {
		allShaders = allShaders.concat(evInvShaders);
	});
	allShaders = allShaders.concat(roeComputeDeltaQTildeShader);
	allShaders = allShaders.concat(roeComputeFluxSlopeShader);
	$.each(roeComputeFluxShader, function(k, fluxShaders) {
		allShaders = allShaders.concat(fluxShaders);
	});
	allShaders.push(roeUpdateStateShader);

	//common to both
	allShaders.push(minReduceShader);

	//display
	$.each(drawToScreenShader, function(k, drawShader) {
		allShaders.push(drawShader);
	});


	//init hydro after gl

	hydro = new Hydro();

	//init controls after hydro
	
	panel = $('#panel');	
	var panelContent = $('#content');
	$('#menu').click(function() {
		if (panelContent.css('display') == 'none') {
			panelContent.show();
		} else {
			panelContent.hide();
		}
	});

	$('#resetSod').click(function(){ hydro.state.resetSod(); });
	$('#resetSodCylinder').click(function(){ hydro.state.resetSodCylinder(); });
	$('#resetWave').click(function(){ hydro.state.resetWave(); });
	$('#resetKelvinHemholtz').click(function(){ hydro.state.resetKelvinHemholtz(); });

	$('#useNoise').change(function() {
		useNoise = $(this).is(':checked');
	});

	(function(){
		var select = $('#gridsize');
		$.each([32, 64, 128, 256, 512, 1024, 2048], function(i,gridsize){
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
	buildSelect('advectMethod', 'advectMethod', advectMethods);
	buildSelect('drawToScreenMethod', 'drawToScreenMethod', drawToScreenMethods);

	$.each([
		'externalForceX',
		'externalForceY'
	], function(i,varName) {
		$('#'+varName).val(window[varName]);
		$('#'+varName).change(function() {
			var v = Number($(this).val());
			if (v !== v) return;	//NaN
			window[varName] = v;
		
			$.each(allShaders, function(k,shader) {
				shader
					.use()
					.setUniforms({
						externalForce : [externalForceX, externalForceY]
					})
					.useNone();
			});
		});
	});
	
	hydro.lastDataMin = Number($('#dataRangeFixedMin').val());
	hydro.lastDataMax = Number($('#dataRangeFixedMax').val());
	hydro.updateLastDataRange = false;
	$('#dataRangeScaleNormalized').change(function() {
		if (!$(this).is(':checked')) return;
		hydro.updateLastDataRange = true;
	});
	$('#dataRangeScaleFixed').change(function() {
		if (!$(this).is(':checked')) return;
		hydro.updateLastDataRange = false;
		hydro.lastDataMin = Number($('#dataRangeFixedMin').val()); 
		hydro.lastDataMax = Number($('#dataRangeFixedMax').val()); 
	});
	$('#dataRangeFixedMin').change(function() {
		if (hydro.updateLastDataRange) return;
		hydro.lastDataMin = Number($('#dataRangeFixedMin').val()); 
	});
	$('#dataRangeFixedMax').change(function() {
		if (hydro.updateLastDataRange) return;
		hydro.lastDataMax = Number($('#dataRangeFixedMax').val()); 
	});
	
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

	//init gl stuff

	glutil.view.ortho = true;
	glutil.view.zNear = -1;
	glutil.view.zFar = 1;
	glutil.view.fovY = 125 / 200 * (xmax - xmin);
	glutil.view.pos[0] = .5;
	glutil.view.pos[1] = .5;

	colorSchemes.Heat = new glutil.GradientTexture({
		width : 256, 
		colors : [
			[0, 0, 0],
			[0, 0, 1],
			[1, 1, 0],
			[1, 0, 0],
		],
		dontRepeat : true
	});

	var isobarSize = 16;
	var isobarData = new Uint8Array(isobarSize);
	for (var i = 1; i < isobarSize; i += 2) {
		isobarData[i] = 255;
	}
	colorSchemes['B&W'] = new glutil.Texture2D({
		width : isobarSize,
		height : 1,
		format : gl.LUMINANCE,
		internalFormat : gl.LUMINANCE,
		data : isobarData,
		minFilter : gl.LINEAR,
		magFilter : gl.NEAREST,
		generateMipmap : true,
		wrap : {s : gl.REPEAT, t : gl.CLAMP_TO_EDGE }
	});

	currentColorScheme = colorSchemes.Heat;

	for (var k in colorSchemes) {
		(function(){
			var v = colorSchemes[k];
			$('#colorScheme').append($('<option>', {
				value : k,
				text : k
			}));
		})();
	}
	$('#colorScheme').change(function() {
		var k = $(this).val();
		currentColorScheme = colorSchemes[k];
	});

	//http://lab.concord.org/experiments/webgl-gpgpu/webgl.html
	for (var channel = 0; channel < 4; ++channel) {
		encodeShader[channel] = new glutil.KernelShader({
			code : mlstr(function(){/*
float shift_right(float v, float amt) {
	v = floor(v) + 0.5;
	return floor(v / exp2(amt));
}

float shift_left(float v, float amt) {
	return floor(v * exp2(amt) + 0.5);
}

float mask_last(float v, float bits) {
	return mod(v, shift_left(1.0, bits));
}

float extract_bits(float num, float from, float to) {
	from = floor(from + 0.5);
	to = floor(to + 0.5);
	return mask_last(shift_right(num, from), to - from);
}

vec4 encode_float(float val) {
	if (val == 0.0)
		return vec4(0, 0, 0, 0);
	float sign = val > 0.0 ? 0.0 : 1.0;
	val = abs(val);
	float exponent = floor(log2(val));
	float biased_exponent = exponent + 127.0;
	float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0;
	
	float t = biased_exponent / 2.0;
	float last_bit_of_biased_exponent = fract(t) * 2.0;
	float remaining_bits_of_biased_exponent = floor(t);
	
	float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0;
	float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0;
	float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0;
	float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0;
	return vec4(byte4, byte3, byte2, byte1);
}

void main() {
	vec4 data = texture2D(tex, pos);
	gl_FragColor = encode_float(data[$channel]);
}
*/}).replace(/\$channel/g, channel),
			texs : ['tex']
		});
	}

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
			glutil.view.pos[0] -= dx / canvas.width * 2 * (aspectRatio * glutil.view.fovY);
			glutil.view.pos[1] += dy / canvas.height * 2 * glutil.view.fovY;
			glutil.updateProjection();
		},
		zoom : function(zoomChange) {
			dragging = true;
			var scale = Math.exp(-zoomFactor * zoomChange);
			glutil.view.fovY *= scale 
			glutil.updateProjection();
		}
	});
	
	//start it off
	$(window).resize(onresize);
	onresize();
	update();


	//check ...	
	var checkCoordAccuracyShader = new glutil.KernelShader({
		code : mlstr(function(){/*
void main() {
	gl_FragColor = vec4(gl_FragCoord.xy*dpos, dpos);
}
*/}),
		uniforms : {
			dpos : ['vec2', [
				1/hydro.state.nx,
				1/hydro.state.nx
			]]
		}
	});

	var nx = hydro.state.nx;
	var tmpTex = new FloatTexture2D(nx, nx); 
	fbo.setColorAttachmentTex2D(0, tmpTex);
	fbo.draw({
		callback : function() {
			gl.viewport(0, 0, tmpTex.width, tmpTex.height);
			quadObj.draw({
				shader : checkCoordAccuracyShader
			});
		}
	});
	var coords = [];
	var maxErr = [];
	for (var channel = 0; channel < 4; ++channel) {
		coords[channel] = getFloatTexData({srcTex:tmpTex, channel:channel});
		maxErr[channel] = 0;
	}
	for (var j = 0; j < nx; ++j) {
		for (var i = 0; i < nx; ++i) {
			var target = [
				(i + .5) / nx,
				(j + .5) / nx,
				1/nx,
				1/nx
			];
			for (var channel = 0; channel < 4; ++channel) {
				var err = Math.abs(coords[channel][i+nx*j] - target[channel]);
				maxErr[channel] = Math.max(err, maxErr[channel]);
			}
		}
	}
	console.log('max coord errors:', maxErr);

});
