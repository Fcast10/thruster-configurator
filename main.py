from fastapi import FastAPI
from pydantic import BaseModel
import math

app = FastAPI()

class InputData(BaseModel):
    dry_mass: float
    delta_v: float
    isp: float = 2300
    rho_g: float = 1.2

class OutputData(BaseModel):
    propellant_mass: float
    wet_mass: float
    tank_vol: float

@app.post("/calculate", response_model=OutputData)
def calculate(data: InputData):
    g0 = 9.80665

    mass_ratio = math.exp(data.delta_v / (data.isp * g0))
    prop_mass = data.dry_mass * (mass_ratio - 1)
    wet_mass = data.dry_mass + prop_mass
    vol = prop_mass / data.rho_g

    return OutputData(
        propellant_mass=prop_mass,
        wet_mass=wet_mass,
        tank_vol=vol
    )

@app.get("/")
def root():
	return{"status": "OK", "message": "Thruster configurator backend is running"}