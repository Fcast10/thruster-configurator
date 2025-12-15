from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import math

app = FastAPI()

app.add_middleware( 
    CORSMiddleware, 
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"], 
)


class InputData(BaseModel):
    dry_mass: float         # [kg]
    delta_v: float          # [m/s]

class OutputData(BaseModel):
    propellant_mass: float  # [kg]
    total_mass: float         # [kg]
    height_cyl: float       # [mm]

# Endpoint di test #

@app.get("/")
def root():
	return{"message":"Thruster configurator API is running"}

# Calcolo principale #

@app.post("/calculate", response_model=OutputData)
def calculate(data: InputData):

    # Fixed Data

    g0 = 9.80665                  # [m/s^2]
    isp = 180                     # [s]
    rho_g = 1000                  # [kg/m^3]
    diameter = 94.2               # [mm]
    radius = diameter/2           # [mm]
    H_dome = 17                   # [mm]
    total_dry_prop_mass = 1.457   # [kg]
    rho_al = 2810                 # [kg/m^3]
    thick_tank = 1                # [mm], better estimate(?)
 

    mass_ratio = math.exp(data.delta_v / (isp * g0))   # dimensionless
    m_f = data.dry_mass + total_dry_prop_mass          # [kg]
    prop_mass = m_f * (mass_ratio - 1)                 # [kg]
    total_vol = prop_mass / rho_g                      # [m^3]
    total_vol = total_vol*1e9                          # [mm^3]
    cross_area = math.pi*(radius)**2                   # [mm^2]
    dome_vol = 2/3*math.pi*(radius)**2*H_dome          # [mm^3]
    cyl_vol = total_vol - 2*dome_vol                   # [mm^3]

    if cyl_vol < 0:
         H_cylinder = 0
    else:
        H_cylinder = cyl_vol/cross_area    # [mm]

    # Preliminary tank sizing, approx: tank geometry, thin walls

    e = math.sqrt(1-(H_dome**2/radius**2))
    A_int = 2*math.pi*radius*H_cylinder + 2*math.pi*radius**2*(1+(1-e**2)/e*math.atanh(e))  # [mm^2]
    tank_mass = rho_al*thick_tank*A_int*1e-9  # [kg]
    tot_mass = m_f + tank_mass + prop_mass 

    return OutputData(
        propellant_mass=prop_mass,
        total_mass=tot_mass,
        height_cyl=H_cylinder
    )


