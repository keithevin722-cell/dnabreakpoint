// Supabase > Project Settings > API. The anon key is public by design; row-level security in supabase/schema.sql protects the data.
// Leave both empty to disable sign-in, visit counting and chat.
const SUPABASE={url:"",anonKey:""};
// Firebase > Project settings > General > Your apps > Web app config. These values are public by design; firestore.rules protects the data.
// Chat stays off until apiKey and appId are filled in.
const FIREBASE={apiKey:"",authDomain:"dna-breakpoint.firebaseapp.com",projectId:"dna-breakpoint",appId:""};
