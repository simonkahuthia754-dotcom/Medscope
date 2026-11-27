import json, traceback
import humanbio
def q(s, v=None):
    try: return humanbio._gql(s, v or {})
    except BaseException: return traceback.format_exc()[-400:]
r = q('query($q:String!){search(queryString:$q,entityNames:["target"],page:{index:0,size:1}){hits{id}}}', {"q": "EGFR"}); print(r)
i = r["data"]["search"]["hits"][0]["id"]
print(json.dumps(q("query($id:String!){target(ensemblId:$id){baselineExpression{rows{tissueBiosampleFromSource tissueBiosampleId median unit datatypeId datasourceId}}}}", {"id": i}))[:1800])
print(json.dumps(q('{__type(name:"BaselineExpression"){fields{name type{name kind ofType{name}}}}}'))[:800])
print(json.dumps(q('{__type(name:"Target"){fields{name}}}'))[:1500])
