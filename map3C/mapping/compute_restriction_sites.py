import re

class ComputeRestrictionSites:

    def has_wildcard(self, pattern):
        # Input pattern can be a list or a string.
        wildcards = re.compile(r'[NMRWYSKHBVD]')
        
        if (isinstance(pattern, list)):
            for p in pattern:
                if re.search(wildcards, p):
                    return True
        else:
            if re.search(wildcards, pattern):
                return True
                
        return False

    def pattern2regexp(self, pattern):
        # Input pattern can be a list or a string.
    
        wildcards = {
            'N': '[ACGT]',
            'M': '[AC]',
            'R': '[AG]',
            'W': '[AT]',
            'Y': '[CT]',
            'S': '[CG]',
            'K': '[GT]',
            'H': '[ACT]',
            'B': '[CGT]',
            'V': '[ACG]',
            'D': '[AGT]',
        }
        
        if isinstance(pattern, list):
            return [ self.pattern2regexp(p) for p in pattern ]
            
        pattern = pattern.upper()
        
        for p, r in wildcards.items():
            pattern = re.sub(p, r, pattern)
            
        return re.compile(pattern.upper())

    def get_match_func(self, pattern):
    
        # Input pattern can be a list or a string.
        
        if not isinstance(pattern, list):
            pattern = [ pattern ]
        
        if self.has_wildcard(pattern):
            pattern = self.pattern2regexp(pattern)
            
            if len(pattern) == 1: # There is only a single pattern.
                pattern = pattern[0] # Use the only element from the list as a single regexp.
                
                def match_single_regexp(segment):
                    if re.match(pattern, segment):
                        return True
                    return False
                    
                return match_single_regexp
            
            else: # There are multiple patterns.
                def match_multi_regexp(segment):
                    for p in pattern:
                        if re.match(p, segment):
                            return True
                    return False
                    
                return match_multi_regexp
            
        else: # No wildcard in any of the patterns.
            
            if len(pattern) == 1: # There is only a single pattern.
                
                pattern = pattern[0] # Use the only element from the list as a single string.
                
                def match_single_string(segment):
                    if segment.startswith(pattern):
                        return True
                    return False
                    
                return match_single_string
            
            else: # There are multiple patterns.
                
                def match_multi_string(segment):
                    for p in pattern:
                        if segment.startswith(p):
                            return True
                    return False
                
                return match_multi_string
            
    def process_input(self):
        f = open(self.reference, 'r')
        g = open(self.output, 'w')
    
        minsize = min([ len(p) for p in self.cut_seqs ])
        maxsize = max([ len(p) for p in self.cut_seqs ])
        matches = self.get_match_func(self.cut_seqs)
    
        segment = ''
        counter = 0
        endl    = ''
        
        for line in f:
    
            line = line.strip()
        
            if line.startswith('>'):
        
                # This is the beginning of a new sequence, but before starting it we must
                # finish processing of the remaining segment of the previous sequence.
        
                while len(segment) > minsize:
                    segment = segment[1:]
                    if matches(segment):
                        g.write(' ' + str(counter - len(segment) + 1))
    
                if counter > 0:
                    g.write(' ' + str(counter)) # Close the previous sequence here.

                firststr=re.split(r'\s+',line[1:])
                g.write(endl+firststr[0])

                segment = ''
                counter = 0
                endl    = '\n'
                
                continue
    
            # Process next line of the sequence.
            
            line = line.upper()
            
            for symbol in line:

                counter += 1
                segment += symbol
                
                while len(segment) > maxsize:
                    segment = segment[1:]

                # Do pattern matching only if segment size equals maxsize.

                if len(segment) == maxsize:
                    if matches(segment):
                        g.write(' ' + str(counter - maxsize + 1)) # maxsize == len(segment)
    
        # Finish the last sequence.
    
        while len(segment) > minsize:
            segment = segment[1:]
            if matches(segment):
                g.write(' ' + str(counter - len(segment) + 1))
    
        if counter > 0:
            g.write(' ' + str(counter))
    
        g.write('\n') # End the output file with a newline.
    
        # Close files.
        
        g.close()
        f.close()        

    def __init__(self, cut_seqs, reference, output):

        self.cut_seqs = [p.upper() for p in cut_seqs]
        self.reference = reference
        self.output = output
        self.process_input()
