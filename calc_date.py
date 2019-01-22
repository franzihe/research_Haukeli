#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import datetime
import calendar


# In[ ]:


def get_dayname(year, mon, day):
    yr = int(year)
    mo = int(mon)
    dy = int(day)
    my_date = datetime.date(yr,mo,dy)
    calday = calendar.day_name[my_date.weekday()]
    calmon = calendar.month_abbr[mo]

    return(calday, calmon);

