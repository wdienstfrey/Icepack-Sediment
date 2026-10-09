"""Check regional input units, flotation, velocity rotation and nodata rejection."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import rasterio
from rasterio.transform import from_origin

ROOT=Path(__file__).resolve().parents[1]

class PreprocessingTests(unittest.TestCase):
    def prepare(self, directory, surface, azimuth=0, missing_bed=False):
        rasters={}
        for name,value in [('bed',-300),('surface',surface),('smb',917),('ux',1000),('uy',0)]:
            path=directory/f'{name}.tif'
            values=np.full((7,7),value,dtype='float64')
            if name=='bed' and missing_bed: values[:]=-9999
            with rasterio.open(path,'w',driver='GTiff',width=7,height=7,count=1,
                               dtype='float64',crs='EPSG:3413',nodata=-9999,
                               transform=from_origin(-3500,3500,1000,1000)) as raster:
                raster.write(values,1)
            rasters[name]={'path':str(path),'scale':1/917 if name=='smb' else 1,
                           'source':'synthetic preprocessing test'}
        manifest={'origin_easting':0,'origin_northing':0,
                  'downstream_azimuth_deg_from_east':azimuth,'maximum_length_m':2000,
                  'width_m':2000,'spacing_m':1000,'rasters':rasters}
        (directory/'manifest.json').write_text(json.dumps(manifest))
        return subprocess.run([sys.executable,str(ROOT/'scripts/prepare_inputs.py'),
            str(directory/'manifest.json'),'--output',str(directory/'inputs.npz')],
            capture_output=True,text=True)

    def test_grounding_units_and_rotation(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            result=self.prepare(directory,100,90)
            self.assertEqual(result.returncode,0,result.stderr)
            with np.load(directory/'inputs.npz') as data:
                np.testing.assert_allclose(data['thickness'],400)
                np.testing.assert_allclose(data['smb'],1)
                np.testing.assert_allclose(data['ux'],0,atol=1e-10)
                np.testing.assert_allclose(data['uy'],-1000)

    def test_floating_thickness(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            result=self.prepare(directory,20)
            self.assertEqual(result.returncode,0,result.stderr)
            with np.load(directory/'inputs.npz') as data:
                np.testing.assert_allclose(data['thickness'],20/(1-917/1024))

    def test_nodata_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            result=self.prepare(Path(folder),100,missing_bed=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('corridor extends beyond valid raster coverage',result.stderr)

if __name__=='__main__': unittest.main()
