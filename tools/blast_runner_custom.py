from tools.blast_runner import BlastRunner



class BlastRunnerCustom(BlastRunner):
    def __init__(self, *args, custom_headers="", **kwargs):
        super(BlastRunnerCustom, self).__init__(*args, **kwargs)
        self.column_headers = self.column_headers + " " + custom_headers

    # def append_column_headers(self, headers_str):
    #     self.custom_headers = self.custom_headers + " " + headers_str

    def get_column_headers(self):
        return self.column_headers