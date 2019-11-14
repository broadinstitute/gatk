"""
To build curl commands from copy pasted forms from the biobank website
"""


FORM_TEXT = """
"""


NAME = "DOWNLOAD.enc"  # Downloaded file's name


test = """
<form name="fetch" action="https://biota.osc.ox.ac.uk/dataset.cgi" method="post">
<input type="hidden" name="id" value="AAA">
<input type="hidden" name="s" value="BBB">
<input type="hidden" name="t" value="CCC">
<input type="hidden" name="i" value="DDD">
<input type="hidden" name="v" value="EEE">
<input class="sub_go" type="submit" value="Fetch">
</form>
"""


def get_fields(txt):
    i = txt.find('''name="fetch"''')
    if i == -1:
        print('Fetch form not in text')
        return
    action, i = get_field(txt, i, '''action="''')
    fields = {'action': action}
    for field in ['id', 's', 't', 'i', 'v']:
        fields[field], i = get_field(txt, i)
    return fields


def get_field(txt, start, target='''value="'''):
    start = txt.find(target, start)
    end = txt.find('''"''', start + len(target))
    return txt[start + len(target): end], end


def fields_to_curl(name, action, id, s, t, i, v):
    return f"""
curl -d "id={id}&s={s}&t={t}&i={i}&v={v}&submit=Fetch" \
-X POST {action} \
-o {name}    
    """


def txt_to_curl(name, txt):
    return fields_to_curl(name, **get_fields(txt))


print(txt_to_curl(NAME, FORM_TEXT))
