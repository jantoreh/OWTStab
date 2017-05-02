function phi = TorsionalModeFun(s,ModeNr)
phi = sin((2*ModeNr-1)*pi*s/(2*s(end)));