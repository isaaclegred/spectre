// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <hdf5.h>
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include "IO/H5/VolumeData.hpp"

#include "ErrorHandling/Assert.hpp"
#include "ErrorHandling/Error.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/Header.hpp"
#include "IO/H5/Helpers.hpp"
#include "IO/H5/Version.hpp"
#include "Utilities/MakeString.hpp"
namespace h5 {
/// We maintain a list of bases and quadratures which are compatible
/// with IO, this allows for efficient storage of Spectral::Basis in
/// volume data files
const std::vector<Spectral::Basis> IO_bases{Spectral::Basis::Chebyshev,
                                            Spectral::Basis::Legendre,
    Spectral::Basis::FiniteDifference};
/// A list of the names of the Bases, for use in a dictionary
const std::vector<std::string> IO_basis_names = []() {
  std::vector<std::string> local_basis_names;
  for (auto& basis : IO_bases) {
    local_basis_names.push_back(MakeString{} << basis);
  }
  return local_basis_names;
}();

/// We maintain a list of quadratures which are compatible
/// with IO, this allows for efficient storage of Spectral::Quadrature in
/// volume data files
const std::vector<Spectral::Quadrature> IO_quadratures{
    Spectral::Quadrature::Gauss, Spectral::Quadrature::GaussLobatto,
    Spectral::Quadrature::CellCentered, Spectral::Quadrature::FaceCentered};
/// A list of the names of the Quadratures, for use in a dictionaray
const std::vector<std::string> IO_quadrature_names = []() {
  std::vector<std::string> local_quadrature_names;
  for (auto& quadrature : IO_quadratures) {
    local_quadrature_names.push_back(MakeString{} << quadrature);
  }
  return local_quadrature_names;
}();
/// Write a dictionary as an attribute to the volume file, can be used
/// to decode integer sequence as values[i] represents the string
/// value encoded which was enconded with integer i
void write_dictionary(const std::string& dict_name,
                      const std::vector<std::string>& values,
                      const detail::OpenGroup& observation_group) {
  h5::write_to_attribute<std::string>(observation_group.id(), dict_name,
                                      values);
}
/// A dictionary `dict_name` is used to decode the integer vector `decodable`
/// into an vector of strings.  The `dict_name` should correspond to
/// a dictionary written with `h5::write_dictonary`
std::vector<std::string> decode_with_dictionary_name(
    const std::string& dict_name, const std::vector<int>& decodable,
    const detail::OpenGroup& observation_group) {
  const auto dict =
      h5::read_rank1_attribute<std::string>(observation_group.id(), dict_name);

  std::vector<std::string> decoded(decodable.size());
  for (size_t i = 0; i < decodable.size(); i++) {
    decoded[i] = dict[static_cast<size_t>(decodable[i])];
  }
  return decoded;
}

}  // namespace h5
