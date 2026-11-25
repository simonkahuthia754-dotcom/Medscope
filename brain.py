"""brain.py - pluggable AI reasoning for MedScope. Stdlib only.
Backends tried in order (first that works):
  1. Ollama  - local, free, truly unlimited (needs a PC/server running Ollama; set OLLAMA_HOST, MEDSCOPE_MODEL)
  2. Gemini  - free tier, rate-limited (set GEMINI_API_KEY)
  3. Groq    - free tier, rate-limited (set GROQ_API_KEY)
No backend is truly 'unlimited' except a local one. The brain only summarizes data you pass in; it can be wrong."""
import json, os, urllib.request

SYS = ("You are a careful biomedical research assistant. Use ONLY the evidence provided. Cite items by number. "
       "Never claim a cure is proven. Flag uncertainty, evidence stage (preclinical/Ph1/Ph2/Ph3/approved), and what is unknown.")

def _post(url, body, hdr=None, t=120):
    rq = urllib.request.Request(url, json.dumps(body).encode(), {"Content-Type": "application/json", **(hdr or {})})
    with urllib.request.urlopen(rq, timeout=t) as r: return json.loads(r.read().decode())

def ollama(p):
    h = os.environ.get("OLLAMA_HOST", "http://localhost:11434")
    return _post(h + "/api/generate", {"model": os.environ.get("MEDSCOPE_MODEL", "llama3.2"), "system": SYS, "prompt": p, "stream": False})["response"]

def gemini(p):
    k = os.environ["GEMINI_API_KEY"]; m = os.environ.get("GEMINI_MODEL", "gemini-2.0-flash")
    r = _post(f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={k}",
              {"systemInstruction": {"parts": [{"text": SYS}]}, "contents": [{"parts": [{"text": p}]}]})
    return r["candidates"][0]["content"]["parts"][0]["text"]

def groq(p):
    r = _post("https://api.groq.com/openai/v1/chat/completions", {"model": os.environ.get("GROQ_MODEL", "llama-3.3-70b-versatile"),
              "messages": [{"role": "system", "content": SYS}, {"role": "user", "content": p}]}, {"Authorization": "Bearer " + os.environ["GROQ_API_KEY"]})
    return r["choices"][0]["message"]["content"]

def ask(p):
    errs = []
    for n, f in (("ollama", ollama), ("gemini", gemini), ("groq", groq)):
        try: return n, f(p)
        except Exception as e: errs.append(f"{n}: {type(e).__name__}")
    raise RuntimeError("no brain available (" + ", ".join(errs) + "). Install Ollama or set GEMINI_API_KEY / GROQ_API_KEY.")

def brief(disease, items):
    """summarize numbered evidence items (papers/trials/drugs) into a cautious research brief"""
    ev = "\n".join(f"[{i+1}] {x}" for i, x in enumerate(items))
    return ask(f"Disease: {disease}\nEvidence:\n{ev}\n\nWrite: (1) current treatment landscape, (2) most promising leads and their evidence stage, "
               f"(3) key gaps/risks, (4) next real-world experiments to run. Cite [n]. Be concise.")
