from skyfield.api import Topos, load, EarthSatellite
import matplotlib.pyplot as plt

# Load the satellite TLE data
line1 = '1 25544U 98067A   22071.13413785  .00000223  00000-0  75281-5 0  9994'
line2 = '2 25544  51.6469 324.2679 0003545  94.2234 265.9077 15.48946211282809'
satellite = EarthSatellite(line1,line2, name="ISS (ZARYA)")

# Load the ephemeris data
ts = load.timescale()
t = ts.utc(2022, 3, range(1, 31))  # Time range

# Compute the position of the satellite
satellite_positions = satellite.at(t)

# Extract the coordinates
x = satellite_positions.position.km[0]
y = satellite_positions.position.km[1]
z = satellite_positions.position.km[2]

# Plot the orbit
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.set_title('Orbit of Satellite')
plt.show()
