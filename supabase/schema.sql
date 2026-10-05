-- Run once in Supabase > SQL Editor. Then copy Project URL and anon key
-- (Project Settings > API) into the SUPABASE constant in index.html, and add
-- the site URL under Authentication > URL Configuration (Site URL + Redirect URLs).
-- The anon key is public by design; these row-level-security rules are what protect the data.

create table if not exists public.page_hits (
  id bigint generated always as identity primary key,
  visitor_id text not null check (char_length(visitor_id) between 8 and 64),
  created_at timestamptz not null default now()
);

alter table public.page_hits enable row level security;

revoke all on public.page_hits from anon, authenticated;
grant insert (visitor_id) on public.page_hits to anon, authenticated;

-- Browsers may add a hit but can never read or change rows.
drop policy if exists "anyone can record a hit" on public.page_hits;
create policy "anyone can record a hit" on public.page_hits
  for insert to anon, authenticated
  with check (char_length(visitor_id) between 8 and 64);

-- Only aggregate numbers are exposed publicly.
create or replace function public.site_stats()
returns json
language sql
security definer
set search_path = public
as $$
  select json_build_object(
    'views', count(*),
    'visitors', count(distinct visitor_id)
  ) from public.page_hits;
$$;

revoke all on function public.site_stats() from public;
grant execute on function public.site_stats() to anon, authenticated;

-- Community chat: public read of visible messages, anonymous posting, owner moderation.
create table if not exists public.chat_messages (
  id bigint generated always as identity primary key,
  client_id text not null check (char_length(client_id) between 8 and 64),
  name text not null default 'Anonymous' check (char_length(name) between 1 and 40),
  body text not null check (char_length(body) between 1 and 500),
  hidden boolean not null default false,
  created_at timestamptz not null default now()
);

create index if not exists chat_messages_client_time on public.chat_messages (client_id, created_at desc);
create index if not exists chat_messages_time on public.chat_messages (created_at desc);

alter table public.chat_messages enable row level security;

-- client_id stays private; hide a message by setting hidden = true in the dashboard.
revoke all on public.chat_messages from anon, authenticated;
grant select (id, name, body, created_at) on public.chat_messages to anon, authenticated;
grant insert (client_id, name, body) on public.chat_messages to anon, authenticated;

drop policy if exists "read visible messages" on public.chat_messages;
create policy "read visible messages" on public.chat_messages
  for select to anon, authenticated using (not hidden);

drop policy if exists "post messages" on public.chat_messages;
create policy "post messages" on public.chat_messages
  for insert to anon, authenticated with check (not hidden);

-- One message per client every 5 seconds (a speed bump, not real spam protection).
create or replace function public.chat_rate_limit()
returns trigger
language plpgsql
security definer
set search_path = public
as $$
begin
  if exists (select 1 from public.chat_messages
             where client_id = new.client_id and created_at > now() - interval '5 seconds') then
    raise exception 'Please wait a few seconds before posting again';
  end if;
  return new;
end;
$$;

drop trigger if exists chat_rate_limit on public.chat_messages;
create trigger chat_rate_limit before insert on public.chat_messages
  for each row execute function public.chat_rate_limit();
