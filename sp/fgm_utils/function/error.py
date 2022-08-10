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
    def __init__(self, std, message = "Funky FGM, skip collection!"):
        self.std = std
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} - std:{self.std}"        
