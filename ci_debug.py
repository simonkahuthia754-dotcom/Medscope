import json, traceback
import humanbio
def q(s, v=None):
    try: return humanbio._gql(s, v or {})
    except BaseException: return traceback.format_exc()[-300:]
print(json.dumps(q('{__type(name:"Target"){fields{name args{name type{name kind ofType{name}}}}}}')["data"]["__type"]["fields"][33:36])[:600])
print(json.dumps(q('{__type(name:"BaselineExpressionRow"){fields{name type{name kind ofType{name}}}}}'))[:1200])
i = "ENSG00000146648"
print(json.dumps(q("query($id:String!){target(ensemblId:$id){baselineExpression{count rows{tissueBiosampleFromSource median unit datatypeId datasourceId}}}}", {"id": i}))[:1500])
print(json.dumps(q("query($id:String!){target(ensemblId:$id){baselineExpression{count rows{tissueBiosample{biosampleName} median}}}}", {"id": i}))[:900])
