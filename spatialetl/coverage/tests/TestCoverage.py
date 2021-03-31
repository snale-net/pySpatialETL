from unittest import TestCase

from spatialetl.coverage import Coverage
from spatialetl.coverage.io.netcdf.symphonie.v293 import SYMPHONIEReader


class CoverageTest(TestCase):

    def setUp(self):
        reader = SYMPHONIEReader("resources/cropped_grid.nc", "resources/cropped_example.nc")
        self.coverage = Coverage(reader)

    def test_check_bbox_validity(self):
        self.fail()

    def test_create_mpi_map(self):
        self.fail()

    def test_update_mpi_map(self):
        self.fail()

    def test_read_metadata(self):
        self.fail()

    def test_get_x_size(self):
        self.fail()

    def test_get_y_size(self):
        self.fail()

    def test_is_regular_grid(self):
        self.fail()

    def test_read_axis_x(self):
        self.fail()

    def test_read_axis_y(self):
        self.fail()

    def test_find_point_index(self):
        self.fail()

    def test_read_variable_bathymetry(self):
        self.fail()

    def test_read_variable_topography(self):
        self.fail()

    def test_read_variable_mesh_size(self):
        self.fail()

    def test_read_variable_x_mesh_size(self):
        self.fail()

    def test_read_variable_y_mesh_size(self):
        self.fail()

    def test_read_variable_2d_sea_binary_mask(self):
        self.fail()

    def test_read_variable_ha(self):
        self.fail()
