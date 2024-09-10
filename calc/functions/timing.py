def nice_time(time_in_seconds):
    if time_in_seconds < 1:
        result = f'{time_in_seconds:0.3g}' + ' s'
    elif time_in_seconds < 60:
        result = f'{time_in_seconds:0.2g}' + ' s'
    elif time_in_seconds < 3600:
        result = f'{(time_in_seconds/60):0.2g}' + ' min'
    elif time_in_seconds < 3600*24:
        result = f'{(time_in_seconds/3600):0.1g}' + ' hour'
    else:
        result = f'{(time_in_seconds/(3600*24)):0.1g}' + ' day'
    return result