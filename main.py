from fastapi import FastAPI
from fastapi import HTTPException
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional
from pydantic import BaseModel
import math
import numpy as np

app = FastAPI()

app.add_middleware( 
    CORSMiddleware, 
    allow_origins=[
    "https://ohm.space",  # dominio del sito web
    "https://www.ohm.space",
    "http://localhost",   # per test in locale
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"], 
)

class InputData(BaseModel):
    dry_mass: float                     # [kg]
    delta_v: Optional[float] = None     # [m/s]
    I_tot: Optional[float] = None       # [Ns]
    P_user: float      # [W]
    CS_standard: bool = False

class OutputData(BaseModel):
    propellant_mass: float  # [kg]
    total_mass: float         # [kg]
    height_cyl: float       # [mm]
    number_of_tanks: int 
    
# Endpoint di test #

@app.get("/")
def root():
	return{"message":"Thruster configurator API is running"}
 
# Calcolo principale #

@app.post("/calculate", response_model=OutputData)
def calculate(data: InputData):

    # Input data

    dry_mass = data.dry_mass
    delta_v = data.delta_v
    I_tot = data.I_tot
    P_user = data.P_user
    CS = data.CS_standard

    if delta_v is None and I_tot is None:
        raise HTTPException(
            status_code = 400,
            detail = "Devi fornire almeno deltaV o impulso totale"
        )
    
    if delta_v is not None and I_tot is not None:
        raise HTTPException(
            status_code = 400,
            detail = "Devi fornire solo uno tra deltaV e impulso totale"
         )
            
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
    max_height = 200              # [mm]

    # Linear interpolation Isp w/ Power
    P1, isp1 = 35, 120
    P2, isp2 = 50, 180
    m = (isp2 - isp1)/(P2 - P1)
    q = isp1 - m*P1

    isp_user = m*P_user + q
    isp = isp_user

    if delta_v is not None:
        mass_ratio = math.exp(delta_v / (isp * g0))     # dimensionless
        m_f = dry_mass + total_dry_prop_mass            # [kg]
        prop_mass = m_f * (mass_ratio - 1)              # [kg]
    elif I_tot is not None:
        prop_mass = I_tot/(g0*isp)
        m_f = dry_mass + total_dry_prop_mass 
    
    total_vol = prop_mass / rho_g                       # [m^3]
    total_vol = total_vol*1e9                           # [mm^3]
    cross_area = math.pi*(radius)**2                    # [mm^2]
    dome_vol = 2/3*math.pi*(radius)**2*H_dome           # [mm^3]
    cyl_vol = total_vol - 2*dome_vol                    # [mm^3]

    if cyl_vol < 0:
         H_cylinder = 0
    else:
        H_cylinder = cyl_vol/cross_area                 # [mm]

    tot_height = H_cylinder + 2*H_dome + 2*thick_tank   # [mm]

    # Preliminary tank sizing, approx: tank geometry, thin walls

    e = math.sqrt(1-(H_dome**2/radius**2))
    # A_int = 2*math.pi*radius*H_cylinder + 2*math.pi*radius**2*(1+(1-e**2)/e*math.atanh(e))  # [mm^2]

    # CubeSat constraints
    nm_tanks = 1
    vol_per_tank = total_vol/nm_tanks

    if CS and tot_height > max_height:
        nm_tanks = 2
        vol_per_tank = total_vol/nm_tanks
        cyl_vol = vol_per_tank - 2*dome_vol

        if cyl_vol < 0:
            H_cylinder = 0
        else:
            H_per_cyl = cyl_vol/cross_area
            H_cylinder = H_per_cyl

        tot_height_per_tank = H_cylinder + 2*H_dome + 2*thick_tank   # [mm]

        if tot_height_per_tank > max_height:
            raise HTTPException(status_code=422, detail="La capacità richiesta supera il limite di altezza anche con 2 tank.")

    A_int_per_tank =2*math.pi*radius*H_cylinder + 2*math.pi*radius**2*(1+(1-e**2)/e*math.atanh(e)) # [mm^2]
    one_tank_mass = rho_al*thick_tank*A_int_per_tank*1e-9  # [kg]
    tank_mass = one_tank_mass*nm_tanks
    tot_mass = m_f + tank_mass + prop_mass 
    
    return OutputData(
        propellant_mass=prop_mass,
        total_mass=tot_mass,
        height_cyl=H_cylinder,
        number_of_tanks=nm_tanks,
    )