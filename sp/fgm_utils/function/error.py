class CrossTime1Error(Exception):
    """Exception raised for erros in cross time determination

    Attributes:
        stage --
        message --
    """

    def __init__(self, stage, message="Not enough cross times are determined!"):
        self.stage = stage
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"cross time stage {self.stage}: {self.message}"

class cdfError(Exception):
    """Exception raised when reading cdf

     Attributes:
        file -- name of cdf file 
        message -- either cdf not found or cdf not read
    """
    def __init__(self, file, message):
        self.file = file
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} - {self.file}"

class funkyFGMError(Exception):
    """Exception raised when reading cdf

     Attributes:
        std -- std of spin rate 
    """
    def __init__(self, err,std, message = "Funky FGM, skip collection!"):
        self.std = std
        self.err = err
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} - err:{self.err} std:{self.std}"        

class SCreadError(Exception):
    """Exception raised when not sci zone found

    """
    def __init__(self, message = "No sci zone found!"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}" 

class spikeError_spikcrosstime(Exception):
    """Exception raised when the cloest cross zero time after spike is not found

    """
    def __init__(self, ctime_idx_time, message = "No cross zero time found with this spike for"):
        self.ctime_idx_time = ctime_idx_time
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}: {self.ctime_idx_time}!" 

class spikeError_spikespin(Exception):
    """Exception raised when spike t1 t2 is not found
    """
    def __init__(self, ctime_idx_time, message = "spike spins determine error(spike_ctime_idx1/2) for"):
        self.ctime_idx_time = ctime_idx_time
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}: {self.ctime_idx_time}!" 


class spikeError_t1t2(Exception):
    """Exception raised when spike t1 t2 is not found after spike - avg

    """
    def __init__(self, ctime_idx_time, message = "No spike t1/t2 found in avg - spike for"):
        self.ctime_idx_time = ctime_idx_time
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}: {self.ctime_idx_time}!" 



