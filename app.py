import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import pandas as pd
import astropy.coordinates
from math import *
from flask import Flask, Response

app = Flask(__name__)

# Read the Exoplanet list, skip comments
Exoplanetlist = pd.read_csv("NASA Exoplanet Archive-all.csv", skiprows=291)
Exoplanetlist = Exoplanetlist[
    ["pl_name", "gaia_id", "hostname", "glat", "glon", "ra", "dec", "sy_dist"]
]

Exoplanetlist.drop_duplicates(["gaia_id"], inplace=True, keep="last")
Exoplanetlist.set_index("gaia_id", drop=True, inplace=True)


@app.route("/get_stars/<planetname>")
def get_stars(planetname: str):
    hostname = Exoplanetlist[Exoplanetlist["pl_name"] == planetname].hostname[0]
    planetra = Exoplanetlist[Exoplanetlist["pl_name"] == planetname].ra[0]
    planetdec = Exoplanetlist[Exoplanetlist["pl_name"] == planetname].dec[0]
    planetdist = Exoplanetlist[Exoplanetlist["pl_name"] == planetname].sy_dist[0]
    planetglon = Exoplanetlist[Exoplanetlist["pl_name"] == planetname].glon[0]
    planetglat = Exoplanetlist[Exoplanetlist["pl_name"] == planetname].glat[0]
    job = Gaia.launch_job(
        f"""SELECT TOP 5000 gaia_source.designation,gaia_source.source_id,gaia_source.pseudocolour,gaia_source.l,gaia_source.b,gaia_source.grvs_mag,gaia_source.ra,gaia_source.dec,gaia_source.distance_gspphot

        FROM gaiadr3.gaia_source

        WHERE

        CONTAINS(

            POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),

            CIRCLE('ICRS',{planetra},{planetdec},{u.Quantity(degrees(atan(20/planetdist)), u.deg).value})

        )=1  AND  (gaiadr3.gaia_source.distance_gspphot>={planetdist-20} AND gaiadr3.gaia_source.distance_gspphot<={planetdist+20})"""
    )

    df = job.get_results().to_pandas()
    df.rename(
        columns={"distance_gspphot": "distance", "l": "glon", "b": "glat"}, inplace=True
    )
    planetcoordinates = astropy.coordinates.spherical_to_cartesian(
        planetdist, radians(planetglat), radians(planetglon)
    )
    exoplanetdata = pd.DataFrame(
        columns=["source_id", "l", "d", "pseudocolour", "brightness"]
    )

    for i in df.index:
        source_id = df.loc[i, "source_id"]
        pseudocolour = df.loc[i, "pseudocolour"]

        d1 = df.loc[i, "distance"]

        starcoordinates = astropy.coordinates.spherical_to_cartesian(
            d1, radians(df.loc[i, "glat"]), radians(df.loc[i, "glon"])
        )

        difference = [0, 0, 0]

        for i in range(3):

            difference[i] = starcoordinates[i] - planetcoordinates[i]

        perspectivecoordinates = astropy.coordinates.cartesian_to_spherical(
            difference[0], difference[1], difference[2]
        )

        exoplanetbrightness = (
            df.loc[i, "grvs_mag"] / (d1**2) * (perspectivecoordinates[0].value ** 2)
        )

        a = pd.DataFrame(
            {
                "source_id": source_id,
                "mag": perspectivecoordinates[0].value,
                "l": perspectivecoordinates[1].value,
                "d": perspectivecoordinates[2].value,
                "pseudocolour": [pseudocolour],
                "brightness": exoplanetbrightness,
            }
        )

        exoplanetdata = pd.concat([exoplanetdata, a])

    exoplanetdata.set_index("source_id", inplace=True, drop=True)

    exoplanetdata = pd.DataFrame(columns=["source_id", "brightness", "pseudocolour"])
    for i in df.index:
        source_id = df.loc[i, "source_id"]
        pseudocolour = df.loc[i, "pseudocolour"]
        d1 = df.loc[i, "distance"]

        starcoordinates = astropy.coordinates.spherical_to_cartesian(
            d1, radians(df.loc[i, "glat"]), radians(df.loc[i, "glon"])
        )

        difference = [0, 0, 0]
        for i in range(3):
            difference[i] = starcoordinates[i] - planetcoordinates[i]

        exoplanetbrightness = (
            df.loc[i, "grvs_mag"]
            / (d1**2)
            * (
                difference[0].value ** 2
                + difference[1].value ** 2
                + difference[2].value ** 2
            )
        )
        a = pd.DataFrame(
            {
                "source_id": source_id,
                "x": difference[0].value,
                "y": difference[1].value,
                "z": difference[2].value,
                "pseudocolour": [pseudocolour],
                "brightness": exoplanetbrightness,
            }
        )
        exoplanetdata = pd.concat([exoplanetdata, a])
    exoplanetdata.set_index("source_id", inplace=True, drop=True)

    return Response(exoplanetdata.to_csv(header=""), mimetype="text/csv")


if __name__ == "__main__":

    app.run()
