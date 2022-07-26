import numpy as np


def write_forcing_cont(a, outfile="forcing_cont.dat"):
    """
    Writes the file forcing_cont.dat that can be used to specify the form of the continuous forcing in forcing.f90

    Parameters
    ----------
    a : numpy array specifying the continuous forcing. Shape is expected to be (3,nx,ny,nz). Note the order of the spatial indices.
    
    outfile : file into which the array should be written

    Example usage
    -------------
    >>> import pencil as pc
    >>> import numpy as np
    >>> dim = pc.read.dim()
    >>> a = np.ones((3,dim.nx,dim.ny,dim.nz))
    >>> pc.util.write_forcing_cont(a)
    """
    a_ = np.reshape(a, np.shape(a), order="F")

    with open(outfile, "w") as f:
        """
        Documention for list-directed IO: https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc5/index.html

        Apparently tabs cause problems in some Fortran compilers, so we have to use spaces in the output.
        """
        nrec = 6  # number of records per line. Just an arbitrary number.
        for elem, i in zip(np.nditer(a_, order="F"), range(0, np.size(a_))):
            if i != 0 and i % nrec == 0:
                f.write("\n")
            f.write("    {}".format(elem))
        f.write("\n")
