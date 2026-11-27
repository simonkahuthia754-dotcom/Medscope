import json, traceback
import humanbio
def q(s, v=None):
    try: return humanbio._gql(s, v or {})
    except BaseException: return traceback.format_exc()[-300:]
r = q('{__schema{types{name}}}'); print([t["name"] for t in r["data"]["__schema"]["types"] if "xpress" in t["name"] or "Baseline" in t["name"]])
print(json.dumps(q('{__type(name:"BaselineExpression"){fields{name type{kind ofType{kind ofType{kind ofType{name}}}}}}}'))[:600])
for n in [t["name"] for t in r["data"]["__schema"]["types"] if "xpress" in t["name"] or "Baseline" in t["name"]]:
    print(n, json.dumps(q('query($n:String!){__type(name:$n){fields{name}}}', {"n": n}))[:700])
