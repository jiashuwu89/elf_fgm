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
        return f"❌ cross time stage {self.stage}: {self.message}"

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
        return f"❌ {self.message} - {self.file}"

class funkyFGMError(Exception):
    """Exception raised when reading cdf

     Attributes:
        std -- std of spin rate 
    """
    def __init__(self, err,std, message = "[PREPROCESS] Funky FGM, skip collection!"):
        self.std = std
        self.err = err
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"❌ {self.message} - median:{self.err} std:{self.std}"   

class funkyFGMError_len(Exception):
    """Exception when too few data
    """
    def __init__(self, message = "[PREPROCESS] Too few data, skip collection!"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"❌ {self.message}"        

class SCreadError(Exception):
    """Exception raised when not sci zone found

    """
    def __init__(self, message = "No sci zone found!"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"❌ {self.message}" 

class spikeError80_spikcrosstime(Exception):
    """Exception raised when the cloest cross zero time after spike is not found

    """
    def __init__(self, ctime_idx_time, message = "1/80s spike: No cross zero time found with this spike for"):
        self.ctime_idx_time = ctime_idx_time
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"⚠️ {self.message}: {self.ctime_idx_time}!" 


class spikeError80_spikespin(Exception):
    """Exception raised when spike t1 t2 is not found
    t1 is the actually start of the spike, usually happen earlier than the actual ctime spike
    """
    def __init__(self, ctime_idx_time, message = "1/80s spike: Spike spins determine error(spike_ctime_idx1/2) for"):
        self.ctime_idx_time = ctime_idx_time
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"⚠️ {self.message}: {self.ctime_idx_time}!" 


class spikeError80_t1t2(Exception):
    """Exception raised when spike t1 t2 is not found after spike - avg

    """
    def __init__(self, ctime_idx_time, message = "1/80s spike: No spike t1/t2 found in avg - spike for"):
        self.ctime_idx_time = ctime_idx_time
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"⚠️ {self.message}: {self.ctime_idx_time}!" 

class spikeError25_spikcrosstime(Exception):
    """Exception raised when the cloest cross zero time after spike is not found

    """
    def __init__(self, ctime_idx_time, message = "2.5 spike: No cross zero time found with this spike for"):
        self.ctime_idx_time = ctime_idx_time
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}: {self.ctime_idx_time}!" 

class fsp_spike_del_error(Exception):
    """Exception raised when the cloest cross zero time after spike is not found

    """
    def __init__(self, message = "[POSTPROCESS] Fsp spike delete error, too few points left!"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"❌ {self.message}" 


class preproc_resample_error(Exception):
    """Exception raised when preprocessing fail
    """
    def __init__(self, message = "[PREPROCESS] Resample fail!"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"❌ {self.message}" 


class preproc_download_error(Exception):
    """Exception raised when preprocessing fail
    """
    def __init__(self, message = "[PREPROCESS] Download CDF fail!"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"❌ {self.message}" 
    

class postproc_fgs_igrf(Exception):
    """Exception raised when fgs_igrf fails in determining the fsp data
    """
    def __init__(self, message = "[POSTPREPROCESS] fgs_igrf fail!"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"❌ {self.message}" 