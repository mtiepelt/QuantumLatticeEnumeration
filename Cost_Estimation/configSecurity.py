"""
    Expected Security of AES under MAXDEPTH restrictions following
    Jaques, S., Naehrig, M., Roetteler, M., Virdia, F.: Implementing grover oracles for quantum key search on AES and LowMC. pp. 280–310 (2020). https://doi.org/10.1007/978-3-030-45724-2_10
    Table 10, 12
"""
aes_expected_security = {
        381:
            {
                -1: 83,
                9999: 83,
                40: 157 - 40,
                64: 157 - 64,
                96: 83, # NOT 157 - 96 = 61
            },
        406:
            {
                -1: 83,
                9999: 83,
                40: 157 - 40,
                64: 157 - 64,
                96: 83, # NOT 157 - 96 = 61,
            },
        623:
            {
                -1: 115,
                9999: 115,
                40: 221 - 40,
                64: 221 - 64,
                96: 221 - 96,
            },
        873:
            {
                -1: 148,
                9999: 148,
                40: 285 - 40,
                64: 285 - 64,
                96: 285 - 96,
            },
    }

"""
    Expected Security of Kyber based on NIST levels according to
    NIST: Submission requirements and evaluation criteria for the post-quantum cryptography standardization process (2016),
    https://csrc.nist.gov/CSRC/media/Projects/Post-Quantum-Cryptography/documents/call-for-proposals-final-dec-2016.pdf
"""
kyber_expected_security = {
        381: 128,
        406: 128,
        623: 192,
        873: 256
    }


# Kyber Parameters and BKZ blocksizes from
# Aono, Y., Nguyen, P.Q., Seito, T., Shikata, J.: Lower bounds on latticeenu meration with extreme pruning. pp. 608–637 (2018). https://doi.org/10.1007/978-3-319-96881-0_21
# Eq. (16)
Kyber = {
        'kyber512' :
            {
                'n' : 512,
                'q' : 3329,
                'beta' : 406
            },
        'kyber768' :
            {
                'n' : 768,
                'q' : 3329,
                'beta' : 623,
            },
        'kyber1024' :
            {
                'n' : 1024,
                'q' : 3329,
                'beta' : 873
            }
    }