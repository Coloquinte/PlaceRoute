{ lib, stdenv, fetchFromGitHub, cmake, lemon-graph, eigen, boost
}:

stdenv.mkDerivation rec {
  pname = "coloquinte";
  version = "0.3.1";

  meta = with lib; {
    description = "Placement library for electronic circuits";
    homepage    = "https://github.com/Coloquinte/PlaceRoute";
    license     = licenses.mit;
    platforms   = platforms.linux;
  };

  src = fetchFromGitHub {
    owner = "coloquinte";
    repo = "PlaceRoute";
    rev = "${version}";
    hash = "sha256-bPDXaNZCNBM0qiu+46cL/zH/41lwqHPqfqTzJaERgVQ=";
  };

  nativeBuildInputs = [
    cmake
  ];

  buildInputs = [
    lemon-graph
    eigen
    boost
  ];
}
