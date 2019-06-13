#!/bin/env python

# Distributed under the MIT License.
# See LICENSE.txt for details.

import glob
from spectre import DataStructures
from spectre.IO import H5 as spectre_h5


def generate_xdmf(file_prefix, output_filename, start_time, stop_time, stride):
    """
    Generate one XDMF file that ParaView and VisIt can use to load the data
    out of the HDF5 files.
    """
    h5files = [(spectre_h5.H5File(filename, 1), filename)
               for filename in glob.glob(file_prefix + "*.h5")]

    element_data = h5files[0][0].get_vol('/element_data')
    temporal_ids_and_values = \
    [(observation_id, element_data.get_observation_value(
        observation_id)) for observation_id in \
     element_data.list_observation_ids()]
    temporal_ids_and_values.sort(key=lambda x: x[1])

    xdmf_output = "<?xml version=\"1.0\" ?>\n" \
        "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n" \
        "<Xdmf Version=\"2.0\">\n" \
        "<Domain>\n" \
        "<Grid Name=\"Evolution\" GridType=\"Collection\" " \
        "CollectionType=\"Temporal\">\n"

    # counter used to enforce the stride
    stride_counter = 0
    for id_and_value in temporal_ids_and_values:
        if id_and_value[1] < start_time:
            continue
        if id_and_value[1] > stop_time:
            break

        stride_counter += 1
        if (stride_counter - 1) % stride != 0:
            continue

        xdmf_output += "  <Grid Name=\"Grids\" GridType=\"Collection\">\n"
        xdmf_output += "    <Time Value=\"%1.14e\"/>\n" % (id_and_value[1])
        # loop over each h5 file
        # Can't use quite the same structure here in the wrapper class as
        # The wrapper classes don't have objects representing observation_id's,
        # So we'll need to hold on to the .vol file, and get at the grid
        # And the extents one at a time.
        for h5file in h5files:
            current_vol = h5file[0].get_vol('/element_data')
            h5temporal_id = id_and_value[0]
            for grid in current_vol.list_grids(h5temporal_id):
                extents = current_vol.get_extents(h5temporal_id, grid)
                number_of_cells = 1
                for x in extents:
                    number_of_cells *= (x - 1)
                data_item = "        <DataItem Dimensions=\"%d %d %d\" " \
                    "NumberType=\"Double\" Precision=\"8\" Format=\"HDF5\">\n" \
                    % (extents[0], extents[1], extents[2])
                data_item_vec = "        <DataItem Dimensions=\"%d %d %d 3\" " \
                "ItemType = \"Function\" Function = \"JOIN($0,$1,$2)\">\n" \
                    % (extents[0], extents[1], extents[2])
                element_path = "          %s:/element_data.vol/" %(h5file[1]) +\
                "ObservationId%s/%s" % (id_and_value[0], grid)
                xdmf_output += \
                    "    <Grid Name=\"%s\" GrideType=\"Uniform\">\n" % (grid)
                # Write topology information
                xdmf_output += "      <Topology TopologyType=\"Hexahedron\" " \
                    "NumberOfElements=\"%d\">\n" % (number_of_cells)
                xdmf_output += "        <DataItem Dimensions=\"%d %d %d 8\" " \
                    "NumberType=\"Int\" Format=\"HDF5\">\n" % (
                        extents[0] - 1, extents[1] - 1, extents[2] - 1)
                xdmf_output += element_path + "/connectivity\n" \
                    "        </DataItem>\n      </Topology>\n"
                # Write geometry/coordinates
                xdmf_output += "      <Geometry Type=\"X_Y_Z\">\n"
                xdmf_output += data_item + element_path + \
                    "/InertialCoordinates_x\n        </DataItem>\n"
                xdmf_output += data_item + element_path + \
                    "/InertialCoordinates_y\n        </DataItem>\n"
                xdmf_output += data_item + element_path + \
                    "/InertialCoordinates_z\n        </DataItem>\n"
                xdmf_output += "      </Geometry>\n"
                # Write tensor components as scalars
                # Except for vectors, which can be written as vectors.
                components = current_vol.list_tensor_components(
                    h5temporal_id, grid)
                components.remove('InertialCoordinates_x')
                components.remove('InertialCoordinates_y')
                components.remove('InertialCoordinates_z')
                for component in components:
                    # If the component is x component of a vector, write a
                    # vector for all three components
                    if component.endswith("_x"):
                        vector = component[:-2]
                        xdmf_output += "      <Attribute Name=\"%s\" " \
                        "AttributeType=\"Vector\" Center=\"Node\">\n" % (
                            vector)
                        xdmf_output += data_item_vec
                        for index in ["_x", "_y", "_z"]:
                            xdmf_output += data_item + element_path  + \
                            "/%s" %(vector) + index + "\n"   + \
                            "        </DataItem>\n"
                        xdmf_output += "        </DataItem>\n"
                        xdmf_output += "      </Attribute>\n"
                    # If the component is not part of a vector,
                    # write it as a scalar
                    elif(not component.endswith("_y") and \
                         not component.endswith("_z")):
                        xdmf_output += "      <Attribute Name=\"%s\" " \
                        "AttributeType=\"Scalar\" Center=\"Node\">\n" % (
                            component)
                        xdmf_output += data_item + element_path + (
                        "/%s\n" % component) + "        </DataItem>\n"
                        xdmf_output += "      </Attribute>\n"
                    # The component is a y or z component of a vector,
                    # it will be processed with the x component
                    else:
                        continue
                xdmf_output += "    </Grid>\n"


        # close time grid
        xdmf_output += "  </Grid>\n"

    xdmf_output += "</Grid>\n</Domain>\n</Xdmf>"

    for h5file in h5files:
        h5file[0].close()

    with open(output_filename + ".xmf", "w") as xmf_file:
        xmf_file.write(xdmf_output)


def parse_args():
    """
    Parse the command line arguments
    """
    import argparse as ap
    parser = ap.ArgumentParser(
        description="Generate XDMF file for visualizing SpECTRE data. "
        "To load the XDMF file in ParaView you must choose the 'Xdmf Reader', "
        "NOT 'Xdmf3 Reader'",
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--file-prefix',
        required=True,
        help="The common prefix of the H5 volume files to load")
    parser.add_argument(
        '--output',
        required=True,
        help="Output file name, an xmf extension will be added")
    parser.add_argument(
        "--stride",
        default=1,
        type=int,
        help="View only every stride'th time step")
    parser.add_argument(
        "--start-time",
        default=0.0,
        type=float,
        help="The earliest time at which to start visualizing. The start-time "
        "value is included.")
    parser.add_argument(
        "--stop-time",
        default=1e300,
        type=float,
        help="The time at which to stop visualizing. The stop-time value is "
        "not included.")
    return parser.parse_args()


if __name__ == "__main__":
    input_args = parse_args()
    generate_xdmf(input_args.file_prefix, input_args.output,
                  input_args.start_time, input_args.stop_time,
                  input_args.stride)
